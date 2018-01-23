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

#include <array>
#if BITPIT_ENABLE_MPI==1
#   include <mpi.h>
#endif

#include "bitpit_IO.hpp"
#include "stencil.hpp"
#include "reconstructionStencil.hpp"

using namespace bitpit;

/*!
* Subtest 003
*
* Testing computation of reconstruction polynomial in cartesian 2D configuration
* using cell values only
*/
int subtest_001()
{
    int dim =2;
    double h = 1.;
    std::array<double,3> centre = {{0.,0.,0.}};

    { // test costant reconstruction stencil
        log::cout() << "Testing 2D cartesian constant" << std::endl;
        std::vector<ReconstructionStencil::Condition> equations(1);

        equations[0].id = 0;
        equations[0].size = h;
        equations[0].coordinate = {{0.,0.,0.}};
        equations[0].data = ReconstructionStencil::ReconstructionData::CELL_VALUE;
        equations[0].type = ReconstructionStencil::ReconstructionType::CONSTRAINT;

        ReconstructionStencil stencil(centre,0,dim);
        stencil.compute(equations);

        stencil.display( log::cout() );
    }


    { // test second order reconstruction stencil
        log::cout() << "Testing 2D cartesian linear, using a cross-type stencil " << std::endl;
        std::vector<ReconstructionStencil::Condition> equations(5);

        equations[0].id = 0;
        equations[0].size = h;
        equations[0].coordinate = {{0.,0.,0.}};
        equations[0].data = ReconstructionStencil::ReconstructionData::CELL_VALUE;
        equations[0].type = ReconstructionStencil::ReconstructionType::CONSTRAINT;

        // west neighbour
        equations[1].id = 1;
        equations[1].size = h;
        equations[1].coordinate = {{-h,0.,0.}};
        equations[1].data = ReconstructionStencil::ReconstructionData::CELL_VALUE;
        equations[1].type = ReconstructionStencil::ReconstructionType::MINIMIZE;

        // east neighbour
        equations[2].id = 2;
        equations[2].size = h;
        equations[2].coordinate = {{h,0.,0.}};
        equations[2].data = ReconstructionStencil::ReconstructionData::CELL_VALUE;
        equations[2].type = ReconstructionStencil::ReconstructionType::MINIMIZE;

        // south neighbour
        equations[3].id = 3;
        equations[3].size = h;
        equations[3].coordinate = {{0.,-h,0.}};
        equations[3].data = ReconstructionStencil::ReconstructionData::CELL_VALUE;
        equations[3].type = ReconstructionStencil::ReconstructionType::MINIMIZE;

        // north neighbour
        equations[4].id = 4;
        equations[4].size = h;
        equations[4].coordinate = {{0.,h,0.}};
        equations[4].data = ReconstructionStencil::ReconstructionData::CELL_VALUE;
        equations[4].type = ReconstructionStencil::ReconstructionType::MINIMIZE;

        ReconstructionStencil stencil(centre,1,dim);
        stencil.compute(equations);

        stencil.display( log::cout() );
    }


    { // test second order reconstruction stencil
        log::cout() << "Testing 2D cartesian linear, using S, W and SW neighbours" << std::endl;
        std::vector<ReconstructionStencil::Condition> equations(4);

        equations[0].id = 0;
        equations[0].size = h;
        equations[0].coordinate = {{0.,0.,0.}};
        equations[0].data = ReconstructionStencil::ReconstructionData::CELL_VALUE;
        equations[0].type = ReconstructionStencil::ReconstructionType::CONSTRAINT;

        // west neighbour
        equations[1].id = 1;
        equations[1].size = h;
        equations[1].coordinate = {{-h,0.,0.}};
        equations[1].data = ReconstructionStencil::ReconstructionData::CELL_VALUE;
        equations[1].type = ReconstructionStencil::ReconstructionType::MINIMIZE;

        // south neighbour
        equations[2].id = 2;
        equations[2].size = h;
        equations[2].coordinate = {{0.,-h,0.}};
        equations[2].data = ReconstructionStencil::ReconstructionData::CELL_VALUE;
        equations[2].type = ReconstructionStencil::ReconstructionType::MINIMIZE;

        // south west neighbour
        equations[3].id = 3;
        equations[3].size = h;
        equations[3].coordinate = {{-h,-h,0.}};
        equations[3].data = ReconstructionStencil::ReconstructionData::CELL_VALUE;
        equations[3].type = ReconstructionStencil::ReconstructionType::MINIMIZE;

        ReconstructionStencil stencil(centre,1,dim);
        stencil.compute(equations);

        stencil.display( log::cout() );
    }

    { // test third order reconstruction stencil
        log::cout() << "Testing 2D cartesian quadratic, using a block-type stencil " << std::endl;
        std::vector<ReconstructionStencil::Condition> equations(9);

        equations[0].id = 0;
        equations[0].size = h;
        equations[0].coordinate = {{0.,0.,0.}};
        equations[0].data = ReconstructionStencil::ReconstructionData::CELL_VALUE;
        equations[0].type = ReconstructionStencil::ReconstructionType::CONSTRAINT;

        int id(1);
        for( int i=-1; i<2; ++i){
            for( int j=-1; j<2; ++j){
                int k=0;

                if( i==0 && j==0 && k==0 ){
                    continue;
                }

                equations[id].id = id;
                equations[id].size = h;
                equations[id].coordinate = {{h*i, h*j, h*k}};
                equations[id].data = ReconstructionStencil::ReconstructionData::CELL_VALUE;
                equations[id].type = ReconstructionStencil::ReconstructionType::MINIMIZE;

                ++id;

            }
        }

        ReconstructionStencil stencil(centre,2,dim);
        stencil.compute(equations);

        stencil.display( log::cout() );
    }

    return 0;
}

/*!
* Subtest 002
*
* Testing computation of reconstruction polynomials in cartesian 2D configuration
* using cell values and derivatives
*/
int subtest_002()
{
    int dim = 2;
    double h = 1.;
    std::array<double,3> centre = {{0.,0.,0.}};

    { // test costant reconstruction stencil
        log::cout() << "Testing 2D cartesian linear, with cell values and derivatives" << std::endl;
        std::vector<ReconstructionStencil::Condition> equations(4);

        equations[0].id = 0;
        equations[0].size = h;
        equations[0].coordinate = {{0.,0.,0.}};
        equations[0].data = ReconstructionStencil::ReconstructionData::CELL_VALUE;
        equations[0].type = ReconstructionStencil::ReconstructionType::CONSTRAINT;

        equations[1].id = 1;
        equations[1].size = h;
        equations[1].coordinate = {{-h,0.,0.}};
        equations[1].data = ReconstructionStencil::ReconstructionData::CELL_VALUE;
        equations[1].type = ReconstructionStencil::ReconstructionType::MINIMIZE;

        equations[2].id = 2;
        equations[2].size = h;
        equations[2].coordinate = {{0.,-h,0.}};
        equations[2].data = ReconstructionStencil::ReconstructionData::CELL_VALUE;
        equations[2].type = ReconstructionStencil::ReconstructionType::MINIMIZE;

        equations[3].id = 3;
        equations[3].coordinate = {{0.5*h,0.,0.}};
        equations[3].direction = {{1.,0.,0.}};
        equations[3].data = ReconstructionStencil::ReconstructionData::POINT_DERIVATIVE;
        equations[3].type = ReconstructionStencil::ReconstructionType::CONSTRAINT;

        ReconstructionStencil stencil(centre,1,dim);
        stencil.compute(equations);

        stencil.display( log::cout() );
    }


    return 0;
}

/*!
* Subtest 003
*
* Testing computation of reconstruction polynomials in 1D configuration
*/
int subtest_003()
{
    int dim = 1;
    double h = 1.;
    std::array<double,3> centre = {{0.,0.,0.}};

    { // test costant reconstruction stencil
        log::cout() << "Testing 1D stretched mesh, quadratic, with cell values " << std::endl;
        std::vector<ReconstructionStencil::Condition> equations(3);

        equations[0].id = 0;
        equations[0].size = h;
        equations[0].coordinate = {{0.,0.,0.}};
        equations[0].data = ReconstructionStencil::ReconstructionData::POINT_VALUE;
        equations[0].type = ReconstructionStencil::ReconstructionType::CONSTRAINT;

        equations[1].id = 1;
        equations[1].size = 0.5*h;
        equations[1].coordinate = {{-0.75*h,0.,0.}};
        equations[1].data = ReconstructionStencil::ReconstructionData::POINT_VALUE;
        equations[1].type = ReconstructionStencil::ReconstructionType::MINIMIZE;

        equations[2].id = 2;
        equations[2].size = 2*h;
        equations[2].coordinate = {{1.5*h,0.,0.}};
        equations[2].data = ReconstructionStencil::ReconstructionData::POINT_VALUE;
        equations[2].type = ReconstructionStencil::ReconstructionType::MINIMIZE;

        ReconstructionStencil stencil(centre,2,dim);
        stencil.compute(equations);
        stencil.display( log::cout() );

        std::vector<double> realCoeff(3);
        std::vector<double> pointValues(3);
        std::array<double,3> dist;

        realCoeff[0] = 1.1;
        realCoeff[1] = 2.2;
        realCoeff[2] = 3.3;

        dist = equations[0].coordinate - centre;
        pointValues[0] = realCoeff[0] +realCoeff[1]*dist[0] +0.5*realCoeff[2]*dist[0]*dist[0]; 

        dist = equations[1].coordinate - centre;
        pointValues[1] = realCoeff[0] +realCoeff[1]*dist[0] +0.5*realCoeff[2]*dist[0]*dist[0]; 

        dist = equations[2].coordinate - centre;
        pointValues[2] = realCoeff[0] +realCoeff[1]*dist[0] +0.5*realCoeff[2]*dist[0]*dist[0]; 

        std::vector<double> testCoeff = stencil.computeCoefficients(pointValues);

        log::cout() << " imposed  polynomial " << realCoeff << std::endl;
        log::cout() << " computed polynomial " << testCoeff << std::endl;

    }


    return 0;
}

/*!
* Subtest 004
*
* Testing computation of reconstruction polynomials in 1D configuration
*/
int subtest_004()
{
    int dim = 1;
    double h = 1.;
    std::array<double,3> centre = {{0.,0.,0.}};

    { // test linear reconstruction stencil
        log::cout() << "Testing 1D stretched mesh, quadratic, with cell values " << std::endl;
        std::vector<ReconstructionStencil::Condition> equations(3);

        equations[0].id = 0;
        equations[0].size = h;
        equations[0].coordinate = {{0.,0.,0.}};
        equations[0].data = ReconstructionStencil::ReconstructionData::POINT_VALUE;
        equations[0].type = ReconstructionStencil::ReconstructionType::CONSTRAINT;

        equations[1].id = 1;
        equations[1].size = h;
        equations[1].coordinate = {{-h,0.,0.}};
        equations[1].data = ReconstructionStencil::ReconstructionData::POINT_VALUE;
        equations[1].type = ReconstructionStencil::ReconstructionType::MINIMIZE;

        equations[2].id = 2;
        equations[2].size = h;
        equations[2].coordinate = {{h,0.,0.}};
        equations[2].data = ReconstructionStencil::ReconstructionData::POINT_VALUE;
        equations[2].type = ReconstructionStencil::ReconstructionType::MINIMIZE;

        ReconstructionStencil stencil(centre,1,dim);
        stencil.compute(equations);
        stencil.display( log::cout() );

        std::array<double,3> point = {{0.5*h,0.,0.}};

        bitpit::StencilScalar pointValueStencil = stencil.computePointValueStencil(point);
        pointValueStencil.display( log::cout() );


    }


    return 0;
}
/*!
 * Main program.
 */
int main(int argc, char *argv[])
{
#if BITPIT_ENABLE_MPI==1
    MPI_Init(&argc,&argv);
#else
    BITPIT_UNUSED(argc);
    BITPIT_UNUSED(argv);
#endif

    // Initialize the logger
    log::manager().initialize(log::COMBINED);

    // Run the subtests
    log::cout() << "Testing calculation of reconstruction stencils" << std::endl;

    int status;
    try {
        status = subtest_001();
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

    try {
        status = subtest_002();
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

    try {
        status = subtest_003();
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

    try {
        status = subtest_004();
        if (status != 0) {
            return status;
        }
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif
}
