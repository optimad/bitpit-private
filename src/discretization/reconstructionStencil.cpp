/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2017 OPTIMAD engineering Srl
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of bitpit.
 *
 *  bitpit is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  bitpit is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with bitpit. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/

#include <cassert>
#include <limits>

#include "bitpit_private_lapacke.hpp"

#include "bitpit_patchkernel.hpp"
#include "bitpit_LA.hpp"

#include "reconstructionStencil.hpp"

namespace bitpit {

/*!
   Constructor
   \param[in] centre the polynomial will be centred in this point
   \param[in] order the order of the polynomial, [0] constant, [1] linear, [2] quadratic
   \param[in] dim the number of space dimensions
 */ 
ReconstructionStencil::ReconstructionStencil(const std::array<double,3> &centre, int order, int dim) 
    : m_centre(centre), m_order(order), m_dim(dim)
{
    assert(m_order==0 || m_order==1 || m_order==2);
    assert(m_dim==1 || m_dim==2 || m_dim==3);

}

/*!
   Computes the weights to be used in order to calculate the coefficients
   of the polynomial. The coefficients are such that the conditions decoded
   in equations are enforced.
   \param[in] equations the equations that define the polynomial
 */
void  ReconstructionStencil::compute(const std::vector<Condition> &equations)
{

    // store equations' indices
    int nEquations = equations.size();  
    for( const Condition &equation : equations){
        m_pattern.push_back(equation.id);
    }

    // Compute the number of polynomial coefficients
    int nCoeff = getCoefficientCount();

    // Compute the equality-costrained and  least-squares matrices.
    // and store the pivoting between the order of equations and the 
    // linear system, since the system requires that the first rows
    // corrispond to the least-squares and the last rows to the constraints
    int nConstraints(0) ;
    int nLeastSquares(0) ;
    std::vector<std::vector<double>> A;
    std::vector<std::vector<double>> C;
    std::vector<int> orderConstraints;
    std::vector<int> orderLeastSquares;

    for( int i=0; i<nEquations; ++i) {

        const Condition &equation = equations[i];

        std::vector<double> equationCoefficients;

        if(equation.data == ReconstructionData::POINT_VALUE ){
            equationCoefficients = getPointValueEquation(equation.coordinate);

        } else if(equation.data == ReconstructionData::CELL_VALUE){
            equationCoefficients = getCellAverageEquation(equation.coordinate, equation.size);

        } else if(equation.data == ReconstructionData::POINT_DERIVATIVE){
            equationCoefficients = getPointDerivativeEquation(equation.coordinate, equation.direction);
        }

        if(equation.type==ReconstructionType::CONSTRAINT){
            C.push_back(equationCoefficients);
            orderConstraints.push_back(i);
            ++nConstraints;
        } else if( equation.type==ReconstructionType::MINIMIZE){
            A.push_back(equationCoefficients);
            orderLeastSquares.push_back(i);
            ++nLeastSquares;
        }
    }

    std::vector<int> pivot;
    pivot.swap(orderLeastSquares);
    pivot.insert( pivot.end(), orderConstraints.begin(), orderConstraints.end() );

    // The linear-constrained are introduced in the least-squares problem
    // through Lagrange multipliers. The resulting linear system is
    //
    // | A^t A  C^t | |x     | = |A^t b|
    // | C      0   | |lambda|   |d    |
    // with A and C the least-squares and costraints equations,
    // and b and d are their corresponding RHSs, respectively.
    // x are the coefficients of the polynom and lambda the lagrange multipliers.
    //
    // This system is denoted by S:
    // |S| |x     | = |A^t 0| |b|
    //     |lambda|   |0   I| |d|
    // 
    // The matrices S and S^-1 are symmetric and 
    // only the upper triangles are computed
    int n = nCoeff+nConstraints;
    std::vector<double> S(n*n,0);

    for( int i=0; i<nCoeff; ++i){

        for( int j=i; j<nCoeff; ++j){

            // compute A^t A on the fly
            double ATrasA(0.);
            for( int k=0; k<nLeastSquares; ++k){
                ATrasA += A[k][i] *A[k][j];
            }

            int l = linearIndexColMajor(i,j,n,n);
            S[l] = ATrasA;
        }

        for( int j=nCoeff; j<n; ++j){
            int l = linearIndexColMajor(i,j,n,n);
            S[l] = C[j-nCoeff][i];
        }
    }

    // Matrix is being factorized and its inverse is computed
    int info;
    std::vector<int> ipiv(n);

    info = LAPACKE_dsytrf( LAPACK_COL_MAJOR, 'U', n, S.data(), n, ipiv.data() );
    if( info != 0 ) {
        log::cout() << info << std::endl;
        throw std::runtime_error (" Error in LAPACKE_dpotrf void ReconstructionStencil::compute(...) ");
    }
    info = LAPACKE_dsytri( LAPACK_COL_MAJOR, 'U', n, S.data(), n, ipiv.data() ) ;
    if( info != 0 ) {
        log::cout() << info << std::endl;
        throw std::runtime_error (" Error in LAPACKE_dpotri void ReconstructionStencil::compute(...) ");
    }

    // Finally the final solution reads
    // |x     | = S^-1 |A^t 0| |b| = S^-1 T |b|
    // |lambda|        |0   I| |d|          |d|
    //
    // Since we are interested only in x (the polynomial coefficients)
    // only the first nCoeff rows of the matrix S^-1 T are computed.
    // The stencil are stored according the stored pivoting.
    m_coefficientWeights.resize(nCoeff*nEquations);
    for( int i=0; i<nCoeff; ++i ) {

        for( int j=0; j<nEquations; ++j) {

            double value(0);
            for(int k=0; k<n; ++k){

                if( k<nCoeff && j < nLeastSquares){
                    int l = linearIndexColMajorSymmetric(i,k,n,n,'U');
                    value += S[l] *A[j][k];

                } else if( k>=nCoeff && j >= nLeastSquares && k-nCoeff == j-nLeastSquares){
                    int l = linearIndexColMajorSymmetric(i,k,n,n,'U');
                    value += S[l];
                }

            }

            m_coefficientWeights[ linearIndexColMajor(i,pivot[j],nCoeff,nEquations) ] = value;
        }
    }
}

/*!
  Returns the weight to be used in order to calculate the polynomial coefficient
  \param[in] coeff polynomial coefficient
  \return weights that should multiply the given data
  */
std::vector<double> ReconstructionStencil::getWeights( const Coefficient &coeff) const
{
    int nEquations(m_pattern.size());
    std::vector<double> weights(nEquations,0);

    if( m_dim < getCoefficientMinimumDimension(coeff)) {
        return weights;
    }

    if( m_order < getCoefficientMinimumOrder(coeff) ){
        return weights;
    }

    int i=getIndexFromCoefficient(coeff);
    int nCoeff = getCoefficientCount();

    for( int j=0; j<nEquations; ++j){
        weights[j] = m_coefficientWeights[ linearIndexColMajor(i,j,nCoeff,nEquations) ];
    }

    return weights;
}

/*!
  Prints the stencil to a stream
  \param[in] str stream to be used
  */
void ReconstructionStencil::display( std::ostream &str, double tollerance) const
{

    for( int index = 0; index < getCoefficientCount(); ++index){

        str << getCoefficientFromIndex(index) << " ";

        std::vector<double> weights = getWeights( getCoefficientFromIndex(index));
        std::vector<double>::iterator weightItr = weights.begin(); 

        size_t cellCount = m_pattern.size();
        std::vector<long>::const_iterator cellItr = m_pattern.begin(); 

        for(size_t i=0; i<cellCount; ++i){

            if( std::abs(*weightItr) >= tollerance ){
                str << "(" << *cellItr << "," << *weightItr << ") " ;
            }

            ++cellItr;
            ++weightItr;
        }
        str << std::endl;

    }

}

/*!
  Determines the coefficients of the linear equation that describes the reconstruction of a point value
  \param[in] point the coordinates of the point
  \return coefficients of the linear equation 
  */
std::vector<double> ReconstructionStencil::getPointValueEquation(const std::array<double,3> &point) const
{

    int nCoeff = getCoefficientCount();
    std::vector<double> coeff(nCoeff,0.);

    std::array<double,3> dist = point - m_centre;

    if(m_dim==1){

        for(int i=0; i<nCoeff; ++i){

            switch (i){
                case 0: // p_0
                    coeff[i] = 1.;
                    break;

                case 1: // p_x
                    coeff[i] = dist[0];
                    break;

                case 2: // p_xx
                    coeff[i] = 0.5 *dist[0]*dist[0];
                    break;
            }

        }


    } else if(m_dim==2){

        for(int i=0; i<nCoeff; ++i){

            switch (i){
                case 0: // p_0
                    coeff[i] = 1.;
                    break;

                case 1: case 2: // p_x, p_y
                    coeff[i] = dist[i-1];
                    break;

                case 3: case 4: // p_xx, p_yy
                    coeff[i] = 0.5 *pow(dist[i-3],2) ;
                    break;

                case 5: // p_xy
                    coeff[i] = dist[0]*dist[1];
                    break;
            }

        }


    } else if(m_dim==3) {

        for(int i=0; i<nCoeff; ++i){

            switch (i){
                case 0: // p_0
                    coeff[i] = 1.;
                    break;

                case 1: case 2: case 3: // p_x, p_y, p_z
                    coeff[i] = dist[i-1];
                    break;

                case 4: case 5: case 6: // p_xx, p_yy, p_zz
                    coeff[i] = 0.5 *pow(dist[i-4],2) ;
                    break;

                case 7: // p_xy
                    coeff[i] = dist[0]*dist[1];
                    break;

                case 8: // p_xz
                    coeff[i] = dist[0]*dist[2];
                    break;

                case 9: // p_yz
                    coeff[i] = dist[1]*dist[2];
                    break;
            }

        }
    }

    return coeff;
}

/*!
  Determines the coefficients of the linear equation that describes the reconstruction of the derivative in a point
  \param[in] point the coordinates of the point
  \param[in] direction the direction of the derivative
  \return coefficients of the linear equation 
  */
std::vector<double> ReconstructionStencil::getPointDerivativeEquation(const std::array<double,3> &point, const std::array<double,3> &direction) const
{

    int nCoeff = getCoefficientCount();
    std::vector<double> coeff(nCoeff,0.);

    std::array<double,3> dist = point - m_centre;

    if(m_dim==1){

        for(int i=0; i<nCoeff; ++i){

            switch (i){
                case 0: // p_0
                    coeff[i] = 0;
                    break;

                case 1: // p_x
                    coeff[i] = 1;
                    break;

                case 2: // p_xx
                    coeff[i] = dist[0];
                    break;
            }

        }

    } else if(m_dim==2){

        for(int i=0; i<nCoeff; ++i){

            switch (i){
                case 0: // p_0
                    coeff[i] = 0.;
                    break;

                case 1: case 2: // p_x, p_y
                    coeff[i] = direction[i-1];
                    break;

                case 3: case 4: // p_xx, p_yy
                    coeff[i] = dist[i-3] *direction[i-3];
                    break;

                case 5: // p_xy
                    coeff[i] = dist[0]*direction[1] +dist[1]*direction[0];
                    break;
            }

        }


    } else if(m_dim==3) {

        for(int i=0; i<nCoeff; ++i){

            switch (i){
                case 0: // p_0
                    coeff[i] = 0.;
                    break;

                case 1: case 2: case 3: // p_x, p_y, p_z
                    coeff[i] = direction[i-1];
                    break;

                case 4: case 5: case 6: // p_xx, p_yy, p_zz
                    coeff[i] = dist[i-4] *direction[i-4] ;
                    break;

                case 7: // p_xy
                    coeff[i] = dist[0]*direction[1] +dist[1]*direction[0];
                    break;

                case 8: // p_xz
                    coeff[i] = dist[0]*direction[2] +dist[2]*direction[0];
                    break;

                case 9: // p_yz
                    coeff[i] = dist[1]*direction[2] +dist[2]*direction[1];
                    break;
            }

        }
    }

    return coeff;
}

/*!
  Determines the coefficients of the linear equation that describes the reconstruction of a cell average
  The method works only for ElementType::Voxel and ElementType::Pixel
  \param[in] cellCentre coordinates of the cell centre
  \param[in] cellSize the edge length of the cell
  \return coefficients of the linear equation 
  */
std::vector<double> ReconstructionStencil::getCellAverageEquation( const std::array<double,3> &cellCentre, double cellSize) const
{

    int nCoeff = getCoefficientCount();
    std::vector<double> coeff(nCoeff,0.);

    std::array<double,3> dist = cellCentre - m_centre;

    if(m_dim==1){

        for(int i=0; i<nCoeff; ++i){

            switch (i){
                case 0: // p_0
                    coeff[i] = 1;
                    break;

                case 1: // p_x
                    coeff[i] = dist[0];
                    break;

                case 2: // p_xx
                    coeff[i] = 0.5 *(dist[0]*dist[0] +cellSize*cellSize/12.);
                    break;
            }

        }

    } else if(m_dim==2){

        for(int i=0; i<nCoeff; ++i){

            switch (i){
                case 0: // p_0
                    coeff[i] = 1.;
                    break;

                case 1: case 2: // p_x, p_y
                    coeff[i] = dist[i-1];
                    break;

                case 3: case 4: // p_xx, p_yy
                    coeff[i] = 0.5 *(dist[i-3]*dist[i-3] +cellSize*cellSize/12.)  ;
                    break;

                case 5: // p_xy
                    coeff[i] = dist[0]*dist[1];
                    break;
            }

        }


    } else if(m_dim==3) {

        for(int i=0; i<nCoeff; ++i){

            switch (i){
                case 0: // p_0
                    coeff[i] = 1.;
                    break;

                case 1: case 2: case 3: // p_x, p_y, p_z
                    coeff[i] = dist[i-1];
                    break;

                case 4: case 5: case 6: // p_xx, p_yy, p_zz
                    coeff[i] = 0.5 *(dist[i-4]*dist[i-4] +cellSize*cellSize/12.) ;
                    break;

                case 7: // p_xy
                    coeff[i] = dist[0]*dist[1];
                    break;

                case 8: // p_xz
                    coeff[i] = dist[0]*dist[2];
                    break;

                case 9: // p_yz
                    coeff[i] = dist[1]*dist[2];
                    break;
            }

        }
    }

    return coeff;
}

/*!
  Compute the number of coefficients of the polynom
  \return number of coefficients
  */
int ReconstructionStencil::getCoefficientCount() const
{
    // Compute the number of polynomial coefficients
    int nWeights(0);
    for(int i=0; i<=m_order; ++i){
        nWeights += factorial(m_dim-1+i) /factorial(m_dim-1) /factorial(i);
    }

    return nWeights;
}

/*!
  Compute the number of coefficients of the polynom
  \return number of coefficients
  */
int ReconstructionStencil::getIndexFromCoefficient( const ReconstructionStencil::Coefficient &coeff) const
{

    if(m_dim==1){

        switch (coeff){
            case P_0: 
                return 0;
                break;

            case P_X: 
                return 1;
                break;

            case P_XX: 
                return 2;
                break;

        }

    } else if(m_dim==2){

        switch (coeff){
            case P_0: 
                return 0;
                break;

            case P_X: 
                return 1;
                break;

            case P_Y: 
                return 2;
                break;

            case P_XX: 
                return 3;
                break;

            case P_YY: 
                return 4;
                break;

            case P_XY:
                return 5;
                break;
        }

    } else if(m_dim==3){

        switch (coeff){
            case P_0: 
                return 0;
                break;

            case P_X: 
                return 1;
                break;

            case P_Y: 
                return 2;
                break;

            case P_Z: 
                return 3;
                break;

            case P_XX: 
                return 4;
                break;

            case P_YY: 
                return 5;
                break;

            case P_ZZ:
                return 6;
                break;

            case P_XY:
                return 7;
                break;

            case P_XZ:
                return 8;
                break;

            case P_YZ:
                return 9;
                break;
        }
    }
}

/*!
  Transfors the index used for solving the least-squares system
  into the enum desctibing the coefficient 
  \param[in] i column index used for the least-squares problem
  \return Enum desriptor of the coefficient
  */
ReconstructionStencil::Coefficient ReconstructionStencil::getCoefficientFromIndex( int i ) const
{
    if(m_dim==1){

        switch (i){
            case 0: 
                return ReconstructionStencil::Coefficient::P_0;
                break;

            case 1: 
                return ReconstructionStencil::Coefficient::P_X;
                break;

            case 2: 
                return ReconstructionStencil::Coefficient::P_XX;
                break;
        }

    } else if(m_dim==2){

        switch (i){
            case 0: 
                return ReconstructionStencil::Coefficient::P_0;
                break;

            case 1: 
                return ReconstructionStencil::Coefficient::P_X;
                break;

            case 2: 
                return ReconstructionStencil::Coefficient::P_Y;
                break;

            case 3: 
                return ReconstructionStencil::Coefficient::P_XX;
                break;

            case 4: 
                return ReconstructionStencil::Coefficient::P_YY;
                break;

            case 5:
                return ReconstructionStencil::Coefficient::P_XY;
                break;
        }

    } else if(m_dim==3){

        switch (i){
            case 0: 
                return ReconstructionStencil::Coefficient::P_0;
                break;

            case 1: 
                return ReconstructionStencil::Coefficient::P_X;
                break;

            case 2: 
                return ReconstructionStencil::Coefficient::P_Y;
                break;

            case 3: 
                return ReconstructionStencil::Coefficient::P_Z;
                break;

            case 4: 
                return ReconstructionStencil::Coefficient::P_XX;
                break;

            case 5: 
                return ReconstructionStencil::Coefficient::P_YY;
                break;

            case 6:
                return ReconstructionStencil::Coefficient::P_ZZ;
                break;

            case 7:
                return ReconstructionStencil::Coefficient::P_XY;
                break;

            case 8:
                return ReconstructionStencil::Coefficient::P_XZ;
                break;

            case 9:
                return ReconstructionStencil::Coefficient::P_YZ;
                break;
        }
    }
}

int ReconstructionStencil::getCoefficientMinimumOrder( const Coefficient &coeff) const
{
    switch(coeff){
        case P_0:
            return 0;
            break;

        case P_X: case P_Y: case P_Z:
            return 1;
            break;

        case P_XX: case P_YY: case P_ZZ: case P_XY: case P_XZ: case P_YZ:
            return 2;
            break;
    }
}

int ReconstructionStencil::getCoefficientMinimumDimension( const Coefficient &coeff) const
{
    switch(coeff){
        case P_0:
            return 0;
            break;

        case P_X: case P_XX: 
            return 1;
            break;

        case P_Y: case P_YY: case P_XY:
            return 2;
            break;

        case P_Z: case P_ZZ: case P_XZ: case P_YZ:
            return 3;
            break;
    }
}
/*!
  Computes the factorial of an integer
  \param[in] x integer
  \return factorial
  */
int ReconstructionStencil::factorial(int x) const
{
    return (x==0) ? 1 : x*factorial(x-1);
}

/*!
   Computes the linear index of an element of a matrix stored in
   column-major ordering
   \param[in] rowIndex row index
   \param[in] colIndex column index
   \param[in] rows number of rows of the matrix
   \param[in] columns number of columns of the matrix
   \return linear index
 */
int ReconstructionStencil::linearIndexColMajor(int rowIndex, int colIndex, int rows, int columns) const
{
    return rowIndex +colIndex*rows; 
}

/*!
   Computes the linear index of an element of a matrix stored in
   row-major ordering
   \param[in] rowIndex row index
   \param[in] colIndex column index
   \param[in] rows number of rows of the matrix
   \param[in] columns number of columns of the matrix
   \return linear index
 */
int ReconstructionStencil::linearIndexRowMajor(int rowIndex, int colIndex, int rows, int columns) const
{
    return rowIndex*columns +colIndex; 
}

/*!
   Computes the linear index of an element of a symmetric matrix stored in
   column-major ordering when only either the upper or lower triangle is stored.
   E.g. if the indices (rowIndex,colIndex) correspond to element in the lower
   triangle, but uplo indicates that upper triangle is stored, the symmetric
   index within the upper triangle is returned
   \param[in] rowIndex row index
   \param[in] colIndex column index
   \param[in] rows number of rows of the matrix
   \param[in] columns number of columns of the matrix
   \param[in] uplo either 'U' or 'L' for the upper or the lower triangle
   \return linear index
 */
int ReconstructionStencil::linearIndexColMajorSymmetric(int rowIndex, int colIndex, int rows, int columns, char uplo) const
{

    assert( uplo=='L' || uplo=='U');

    if( (uplo=='U' && colIndex < rowIndex) || (uplo=='L' && colIndex > rowIndex) ){
        return linearIndexColMajor( colIndex, rowIndex, rows, columns );
    } else {
        return linearIndexColMajor( rowIndex, colIndex, rows, columns );
    }
}

/*!
   Computes the linear index of an element of a symmetric matrix stored in
   row-major ordering when only either the upper or lower triangle is stored.
   E.g. if the indices (rowIndex,colIndex) correspond to element in the lower
   triangle, but uplo indicates that upper triangle is stored, the symmetric
   index within the upper triangle is returned
   \param[in] rowIndex row index
   \param[in] colIndex column index
   \param[in] rows number of rows of the matrix
   \param[in] columns number of columns of the matrix
   \param[in] uplo either 'U' or 'L' for the upper or the lower triangle
   \return linear index
 */
int ReconstructionStencil::linearIndexRowMajorSymmetric(int rowIndex, int colIndex, int rows, int columns, char uplo) const
{

    assert( uplo=='L' || uplo=='U');

    if( (uplo=='U' && colIndex < rowIndex) || (uplo=='L' && colIndex > rowIndex) ){
        return linearIndexRowMajor( colIndex, rowIndex, rows, columns );
    } else {
        return linearIndexRowMajor( rowIndex, colIndex, rows, columns );
    }
}

/*!
  Display matrix stored in column major ordering to output stream 
  in a nicely formatted form.
  \param[in,out] out output stream
  \param[in] A matrix to be displayed
  \param[in] rows number of rows
  \param[in] columns number of columns
  */
void ReconstructionStencil::displayColMajor( std::ostream &out, double* A, int rows, int columns ) const
{
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            out << A[linearIndexColMajor(i,j,rows,columns)] << "  ";
        }
        out << std::endl;
    }
}

/*!
  Display matrix stored in row major ordering to output stream 
  in a nicely formatted form.
  \param[in,out] out output stream
  \param[in] A matrix to be displayed
  \param[in] rows number of rows
  \param[in] columns number of columns
  */
void ReconstructionStencil::displayRowMajor( std::ostream &out, double* A, int rows, int columns ) const
{
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            out << A[linearIndexRowMajor(i,j,rows,columns)] << "  " ;
        }
        out << std::endl;
    }
}

/*!
  Display matrix stored in column major ordering to output stream 
  in a nicely formatted form.
  \param[in,out] out output stream
  \param[in] A matrix to be displayed
  \param[in] rows number of rows
  \param[in] columns number of columns
  */
void ReconstructionStencil::displayColMajorSymmetric( std::ostream &out, double* A, int rows, int columns, char uplo ) const
{
    assert( uplo=='L' || uplo=='U');

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            out << A[linearIndexColMajorSymmetric(i,j,rows,columns,uplo)] << "  ";
        }
        out << std::endl;
    }
}

/*!
  Display symmteric matrix stored in row major ordering to output stream 
  in a nicely formatted form.
  \param[in,out] out output stream
  \param[in] A matrix to be displayed
  \param[in] rows number of rows
  \param[in] columns number of columns
  */
void ReconstructionStencil::displayRowMajorSymmetric( std::ostream &out, double* A, int rows, int columns, char uplo ) const
{
    assert( uplo=='L' || uplo=='U');

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            out << A[linearIndexRowMajorSymmetric(i,j,rows,columns,uplo)] << "  ";
        }
        out << std::endl;
    }
}

std::ostream& operator<<(std::ostream &out, const ReconstructionStencil::Coefficient &coeff){
    switch (coeff){
        case ReconstructionStencil::Coefficient::P_0: 
            out << "P_0";
            break;

        case ReconstructionStencil::Coefficient::P_X: 
            out << "P_X";
            break;

        case ReconstructionStencil::Coefficient::P_Y: 
            out << "P_Y";
            break;

        case ReconstructionStencil::Coefficient::P_Z: 
            out << "P_Z";
            break;

        case ReconstructionStencil::Coefficient::P_XX: 
            out << "P_XX";
            break;

        case ReconstructionStencil::Coefficient::P_YY: 
            out << "P_YY";
            break;

        case ReconstructionStencil::Coefficient::P_ZZ:
            out << "P_ZZ";
            break;

        case ReconstructionStencil::Coefficient::P_XY:
            out << "P_XY";
            break;

        case ReconstructionStencil::Coefficient::P_XZ:
            out << "P_XZ";
            break;

        case ReconstructionStencil::Coefficient::P_YZ:
            out << "P_YZ";
            break;
    }
}

}
