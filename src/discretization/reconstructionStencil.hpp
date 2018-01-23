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

#ifndef __BTPIT_RECONSTRUCTIONSTENCIL_HPP__
#define __BTPIT_RECONSTRUCTIONSTENCIL_HPP__

#include <array>
#include <vector>


namespace bitpit {


class ReconstructionStencil {

public:

    enum Coefficient{
        P_0 = 0,
        P_X = 1,
        P_Y = 2,
        P_Z = 3,
        P_XX = 4,
        P_YY = 5,
        P_ZZ = 6,
        P_XY = 7,
        P_XZ = 8,
        P_YZ = 9
    };

    enum ReconstructionData{
        POINT_VALUE,
        POINT_DERIVATIVE,
        CELL_VALUE
    };

    enum ReconstructionType{
        CONSTRAINT,
        MINIMIZE
    };

    struct Condition{
        ReconstructionType type;
        ReconstructionData data;
        long id;
        std::array<double,3> coordinate;
        std::array<double,3> direction;
        double size;
    };

    ReconstructionStencil(const std::array<double,3> &centre, int order, int dimensions);

    void compute( const std::vector<Condition> &equations );
    std::vector<double> computeCoefficients( const std::vector<double> &values );

    std::vector<double> getWeights( const Coefficient &coeff ) const;
    void display(std::ostream &out, double tollerance = 1.e-10) const;

private:
    int m_order; 
    int m_dim; 
    std::array<double,3> m_centre;
    std::vector<long> m_pattern;
    std::vector<double> m_coefficientWeights;

    std::vector<double> getPointValueEquation( const std::array<double,3> &point) const ;
    std::vector<double> getPointDerivativeEquation( const std::array<double,3> &point, const std::array<double,3> &direction) const ;
    std::vector<double> getCellAverageEquation( const std::array<double,3> &cellCentre, double cellSize) const;

    int getCoefficientCount() const;
    int getIndexFromCoefficient( const Coefficient &coeff) const;
    Coefficient getCoefficientFromIndex(int i) const;
    int getCoefficientMinimumOrder( const Coefficient &coeff) const;
    int getCoefficientMinimumDimension( const Coefficient &coeff) const;

    int factorial(int x) const;
    int linearIndexColMajor(int rowIndex, int colIndex, int rows, int columns) const;
    int linearIndexRowMajor(int rowIndex, int colIndex, int rows, int columns) const;
    int linearIndexColMajorSymmetric(int rowIndex, int colIndex, int rows, int columns, char uplo) const;
    int linearIndexRowMajorSymmetric(int rowIndex, int colIndex, int rows, int columns, char uplo) const;
    void displayColMajor( std::ostream &out, double* A, int rows, int columns) const;
    void displayRowMajor( std::ostream &out, double* A, int rows, int columns) const;
    void displayColMajorSymmetric( std::ostream &out, double* A, int rows, int columns, char uplo) const;
    void displayRowMajorSymmetric( std::ostream &out, double* A, int rows, int columns, char uplo) const;
};

std::ostream& operator<<(std::ostream &out, const ReconstructionStencil::Coefficient &coefficient);

}

#endif
