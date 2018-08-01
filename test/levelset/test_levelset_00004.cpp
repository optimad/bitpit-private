
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
# include <ctime>
# include <chrono>

#if BITPIT_ENABLE_MPI==1
# include <mpi.h>
#endif

// bitpit
# include "bitpit_CG.hpp"
# include "bitpit_surfunstructured.hpp"
# include "bitpit_voloctree.hpp"
# include "bitpit_levelset.hpp"


/*!
* Subtest 001
*
* Testing levelset refinement.
*/
int subtest_001()
{
    // 2D test case
    uint8_t dimensions(2);

    // used for timing
    std::chrono::time_point<std::chrono::system_clock> start, end;

    // First Input geometry
    bitpit::log::cout()<< " - Loading stl geometry" << std::endl;
    std::unique_ptr<bitpit::SurfUnstructured> STL0( new bitpit::SurfUnstructured (0,1,dimensions) );

    STL0->importDGF("./data/naca0012.dgf");
    STL0->deleteCoincidentVertices() ;
    STL0->buildAdjacencies() ;

    // Transllate geometry to have its centre in origin
    std::array<double,3> meshMin, meshMax;
    STL0->getBoundingBox(meshMin, meshMax) ;
    std::array<double,3> disalignment = -0.5 *(meshMin+meshMax);
    for( bitpit::Vertex &vertex : STL0->getVertices() ){
        vertex.translate(disalignment);
    }
    STL0->updateBoundingBox(true);


    STL0->getVTK().setName("geometry_004") ;
    STL0->getVTK().setCounter(0);
    STL0->write() ;

    bitpit::log::cout()<< "n. vertex: " << STL0->getVertexCount() << std::endl;
    bitpit::log::cout()<< "n. simplex: " << STL0->getCellCount() << std::endl;

    // Create initial octree mesh
    bitpit::log::cout()<< " - Setting mesh" << std::endl;

    STL0->getBoundingBox( meshMin, meshMax ) ;
    bitpit::log::cout() << meshMin << std::endl;
    bitpit::log::cout() << meshMax << std::endl;

    double h(0);
    for( int i=0; i<3; ++i){
        h = std::max( h, meshMax[i]-meshMin[i] ) ;
    }
    h *= 4;

    double dh = h / 128.;
    bitpit::VolOctree mesh(dimensions, meshMin, h, dh );
    mesh.update() ;

    mesh.getVTK().setName("levelset_004") ;
    mesh.getVTK().setCounter() ;

    // Set levelset configuration
    bitpit::LevelSet levelset;
    levelset.setMesh(&mesh) ;

    int id0 = levelset.addObject(std::move(STL0),BITPIT_PI) ;
    bitpit::LevelSetObject &object0 = levelset.getObject(id0);
    object0.setPropagateSign(true);
    object0.enableVTKOutput(bitpit::LevelSetWriteField::DEFAULT);


    // Compute and write level set on initial mesh
    bitpit::log::cout() << " - Compute initial levelset " << std::endl;
    start = std::chrono::system_clock::now();
    object0.compute();
    end = std::chrono::system_clock::now();
    int elapsed_ini = std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

    mesh.write() ;

    // Move levelset and write
    bitpit::log::cout() << " - Displace surface " << std::endl;

    int nSteps = 100;
    double step = 0.9 *h / (double) nSteps;

    std::array<double,3> translation = {{step, step, 0.}};
    std::array<double,3> centre = {{0., 0., 0.}};
    std::array<double,3> axis = {{0., 0., 1.}};
    double angle =  -2. *BITPIT_PI /(double) nSteps;

    int elapsed_dis(0);
    for( int i=0; i<100; ++i){

        bitpit::log::cout() << i << std::endl;

        start = std::chrono::system_clock::now();
        object0.displaceSurface(translation, centre, axis, angle);
        end = std::chrono::system_clock::now();
        centre += translation;

        elapsed_dis += std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();

        mesh.write() ;
    }

    bitpit::log::cout() << "elapsed time initialization " << elapsed_ini << " ms" << std::endl;
    bitpit::log::cout() << "elapsed time displacement " << elapsed_dis << " ms" << std::endl;

    return 0;
};

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
    bitpit::log::cout() << "Testing levelset langrangian movement" << std::endl;

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
