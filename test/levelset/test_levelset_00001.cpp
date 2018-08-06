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


//Standard Template Library
# define _USE_MATH_DEFINES

# include <cmath>
# include <ctime>
# include <chrono>

# include <array>

#if BITPIT_ENABLE_MPI==1
# include <mpi.h>
#endif

// bitpit
# include "bitpit_levelset.hpp"
# include "bitpit_surfunstructured.hpp"
# include "bitpit_volcartesian.hpp"

/*!
* Generate a 2D surface mesh.
*
* \param[in,out] mesh on output it will contain the generated surface mesh
*/
void Generate2DSurfMesh( bitpit::SurfUnstructured &mesh )
{

    const double R = 1.0;
    const long N = 32;
    double theta;
    double dtheta = 2. * BITPIT_PI/((double) N);
    
    std::array<double,3> point;
    std::vector<long> connect(2, bitpit::Element::NULL_ID);

    // Create vertex list
    // Use non-consecutive vertex ids to test if the levelset can handle them.
    const long vertexIdOffset = 101;
    const long vertexIdStride = 2;

    point[2] = 0.0;
    for (long i = 0; i < N; ++i) {
        theta = ((double) i) * dtheta;
        point[0] = R * cos( theta );
        point[1] = R * sin( theta );
        mesh.addVertex(point, vertexIdOffset + vertexIdStride * i);
    }

    // Create simplex list
    // Use non-consecutive cell ids to test if the levelset can handle them.
    const long cellIdOffset = 202;
    const long cellIdStride = 3;

    for (long i = 0; i < N; ++i) {
        connect[0] = vertexIdOffset + vertexIdStride * i;
        connect[1] = vertexIdOffset + vertexIdStride * ((i + 1) % N);
        mesh.addCell(bitpit::ElementType::LINE, true, connect, cellIdOffset + cellIdStride * i);
    }

}

/*!
* Subtest 001
*
* Testing basic features of a 2D levelset.
*/
int subtest_001()
{
    int dimensions(2) ;

    // Input geometry
    bitpit::log::cout() << " - Loading dgf geometry" << std::endl;
    std::unique_ptr<bitpit::SurfUnstructured> STL( new bitpit::SurfUnstructured(0,1,dimensions) );

    Generate2DSurfMesh( *(STL.get()) ) ;

    STL->getVTK().setName("geometry_001") ;
    STL->write() ;

    bitpit::log::cout() << "n. vertex: " << STL->getVertexCount() << std::endl;
    bitpit::log::cout() << "n. simplex: " << STL->getCellCount() << std::endl;

    // create cartesian mesh around geometry 
    bitpit::log::cout() << " - Setting mesh" << std::endl;
    std::array<int,3> nc = {{64, 64, 0}} ;

    std::array<double,3> meshMin, meshMax;
    STL->getBoundingBox( meshMin, meshMax );

    std::array<double,3> delta = meshMax -meshMin;
    meshMin -=  0.1*delta;
    meshMax +=  0.1*delta;

    delta = meshMax -meshMin;

    bitpit::VolCartesian mesh( 1, dimensions, meshMin, delta, nc);
    mesh.update() ;
    mesh.getVTK().setName("levelset_001") ;


    // Configure levelset
    bitpit::LevelSet levelset ;
    levelset.setMesh(&mesh) ;

    // Object 0 from STL
    int id0 = levelset.addObject(std::move(STL),BITPIT_PI) ;
    bitpit::LevelSetObject &object0 = levelset.getObject(id0);
    object0.enableVTKOutput(bitpit::LevelSetWriteField::VALUE);

    // Object 1 is the envelope of all cells within R=0.5
    std::unordered_set<long> mask;
    for( auto & cell : mesh.getCells() ){
        const long &id = cell.getId() ;
        std::array<double,3> center = mesh.evalCellCentroid(id);
        double r = norm2(center);
        if(r<=0.5){
            mask.insert(id);
        }
    }

    int id1 = levelset.addObject(mask);
    bitpit::LevelSetObject &object1 = levelset.getObject(id1);
    object1.enableVTKOutput(bitpit::LevelSetWriteField::VALUE);

    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
    object0.compute();
    object1.compute();
    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();

    int elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
    bitpit::log::cout() << "elapsed time: " << elapsed_seconds << " ms" << std::endl;

    bitpit::log::cout() << " - Exporting data" << std::endl;
    mesh.write() ;

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
	bitpit::log::manager().initialize(bitpit::log::COMBINED);

	// Run the subtests
	bitpit::log::cout() << "Testing basic levelset features" << std::endl;

	int status;
	try {
		status = subtest_001();
		if (status != 0) {
			return status;
		}
	} catch (const std::exception &exception) {
		bitpit::log::cout() << exception.what();
		exit(1);
	}

#if BITPIT_ENABLE_MPI==1
	MPI_Finalize();
#endif
}
