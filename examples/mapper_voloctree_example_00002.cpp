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

/**
 * \example mapper_voloctree_example_00002.cpp
 * 
 * \brief Mesh mapping computing between adapted voloctree meshes.
 * <b>To run</b>: ./mapper_voloctree_example_00002 \n
 */ 

#include <array>
#if BITPIT_ENABLE_MPI==1
#include <mpi.h>
#endif

#include "bitpit_POD.hpp"
#include "bitpit_patchkernel.hpp"
#include "bitpit_voloctree.hpp"

using namespace bitpit;

/**
 * Run the example.
 */ 
void run()
{

    log::cout() << "  >> 2D octree patch" << "\n";

    /**
     * Create the tree
     */
    double x_0 = 10.;
    double y_0 = 20.;
    double z_0 = 30.;
    double l   = 1.5;

    std::unique_ptr<PabloUniform> treePointer = std::unique_ptr<PabloUniform>(new PabloUniform(x_0, y_0, z_0, l, 2));
    PabloUniform &octree = *treePointer;

    std::cout << " Origin : ( " << octree.getX0() << ", " << octree.getY0() << ", " << octree.getZ0() << " )" << std::endl;
    std::cout << " Length : " << octree.getL() << std::endl;

    /** Refine and write the octree*/
    octree.adaptGlobalRefine();
    octree.adaptGlobalRefine();
    octree.adaptGlobalRefine();
    octree.adaptGlobalRefine();

    /**
     * Create the patch from the existing tree
     */
    /** Create the original patch */
    VolOctree *patch_2D_original = new VolOctree(std::move(treePointer), &treePointer);

    patch_2D_original->update();

#if BITPIT_ENABLE_MPI==1
    /** Partition the patch */
    patch_2D_original->partition(true);
#endif

    /**
     * Create the new tree
     */
    std::unique_ptr<PabloUniform> treePointer2 = std::unique_ptr<PabloUniform>(new PabloUniform(x_0, y_0, z_0, l, 2));
    PabloUniform &octree2 = *treePointer2;

    /** Refine and write the octree */
    octree2.adaptGlobalRefine();
    octree2.adaptGlobalRefine();
    octree2.adaptGlobalRefine();
    octree2.adaptGlobalRefine();

    /** Create a new patch */
    VolOctree *patch_2D = new VolOctree(std::move(treePointer2), &treePointer2);
#if BITPIT_ENABLE_MPI==1
    patch_2D->setVTKWriteTarget(PatchKernel::WriteTarget::WRITE_TARGET_CELLS_INTERNAL);

    /** Partition the patch */
    patch_2D->partition(true);
#endif

    srand(1);

    /** Refine the patch */
    for (int k = 0; k < 4; ++k) {
        long nCells = patch_2D->getCellCount();
        log::cout() << std::endl;
        log::cout() << ">> Marking the cells to adapt... " << std::endl;

        for (int i = 0; i < 100; ++i) {
            long cellId = rand() % nCells * 2;
            if (!patch_2D->getCells().exists(cellId)) {
                continue;
            }

            patch_2D->markCellForRefinement(cellId);
            for (auto neighId : patch_2D->findCellNeighs(cellId)) {
                patch_2D->markCellForRefinement(neighId);
            }
        }

        for (int i = 0; i < 50; ++i) {
            long cellId = rand() % nCells * 2;
            if (!patch_2D->getCells().exists(cellId)) {
                continue;
            }

            patch_2D->markCellForCoarsening(cellId);
            for (auto neighId : patch_2D->findCellNeighs(cellId)) {
                patch_2D->markCellForCoarsening(neighId);
            }
        }

        log::cout() << std::endl;
        log::cout() << ">> Initial number of cells... " << nCells << std::endl;

        patch_2D->update();

        nCells = patch_2D->getCellCount();
        log::cout() << ">> Final number of cells... " << nCells << std::endl;
    }

    /** Show patch info */
    log::cout() << "Cell count:   " << patch_2D->getCellCount() << std::endl;
    log::cout() << "Vertex count: " << patch_2D->getVertexCount() << std::endl;

    patch_2D->getVTK().setName("mesh_random.0");
#if BITPIT_ENABLE_MPI==1
    patch_2D->setVTKWriteTarget(PatchKernel::WriteTarget::WRITE_TARGET_CELLS_INTERNAL);
#endif
    patch_2D->write();

    /** Create mapper object */
    MapperVolOctree mapobject(patch_2D_original, patch_2D);

    /** Map the two meshes */
    mapobject.initialize(true);

    patch_2D_original->getVTK().setName("mesh_original.0");
#if BITPIT_ENABLE_MPI==1
    patch_2D_original->setVTKWriteTarget(PatchKernel::WriteTarget::WRITE_TARGET_CELLS_INTERNAL);
#endif
    patch_2D_original->write();

    std::vector<uint64_t> refineList;
    refineList.push_back(  7);
    refineList.push_back( 13);
    refineList.push_back( 15);
    refineList.push_back( 26);
    refineList.push_back( 27);
    refineList.push_back( 31);
    refineList.push_back( 37);
    refineList.push_back( 39);
    refineList.push_back( 49);
    refineList.push_back( 50);
    refineList.push_back( 51);
    refineList.push_back( 53);
    refineList.push_back( 55);
    refineList.push_back( 63);
    refineList.push_back( 78);
    refineList.push_back(100);
    refineList.push_back(102);
    refineList.push_back(105);
    refineList.push_back(108);
    refineList.push_back(109);
    refineList.push_back(110);
    refineList.push_back(135);
    refineList.push_back(141);
    refineList.push_back(143);
    refineList.push_back(146);
    refineList.push_back(147);
    refineList.push_back(151);
    refineList.push_back(153);
    refineList.push_back(154);
    refineList.push_back(155);
    refineList.push_back(157);
    refineList.push_back(159);
    refineList.push_back(165);
    refineList.push_back(167);
    refineList.push_back(183);
    refineList.push_back(198);
    refineList.push_back(204);
    refineList.push_back(206);
    refineList.push_back(225);
    refineList.push_back(228);
    refineList.push_back(229);
    refineList.push_back(230);

    for (uint64_t ind : refineList) {
#if BITPIT_ENABLE_MPI==1
        int owner = patch_2D_original->getTree().getOwnerRank(ind);
        if (patch_2D_original->getRank() == owner){
            uint32_t lind = patch_2D_original->getTree().getLocalIdx(ind, owner);
#else
            uint32_t lind = patch_2D_original->getTree().getLocalIdx(ind);
#endif
            VolOctree::OctantInfo octinfo(lind, true);
            long id = patch_2D_original->getOctantId(octinfo);
            patch_2D_original->markCellForRefinement(id);
#if BITPIT_ENABLE_MPI==1
        }
#endif
    }

    std::vector<adaption::Info> infoAdapt = patch_2D_original->adaptionPrepare(true);

    mapobject.adaptionPrepare(infoAdapt, true);

    infoAdapt = patch_2D_original->adaptionAlter(true);

    mapobject.adaptionAlter(infoAdapt, true, true);

    patch_2D_original->adaptionCleanup();

    mapobject.adaptionCleanup();

    /** Show patch info */
    log::cout() << "Cell count:   " << patch_2D_original->getCellCount() << std::endl;
    log::cout() << "Vertex count: " << patch_2D_original->getVertexCount() << std::endl;

    /** Define data on original mesh and write */
    PiercedStorage<double> data(1, &patch_2D_original->getCells());
    std::vector<double> vdata(patch_2D_original->getInternalCount());
    int count = 0;
    for (Cell & cell : patch_2D_original->getCells()){
        if (cell.isInterior()){
            long id = cell.getId();
            data[id] = double(patch_2D_original->getCellLevel(id));
            vdata[count] = data[id];
            count++;
        }
    }

    patch_2D_original->getVTK().setName("mesh_original.0");
    patch_2D_original->getVTK().addData("data", VTKFieldType::SCALAR, VTKLocation::CELL, vdata);
#if BITPIT_ENABLE_MPI==1
    patch_2D_original->setVTKWriteTarget(PatchKernel::WriteTarget::WRITE_TARGET_CELLS_INTERNAL);
#endif
    patch_2D_original->write();

    /** Map data on second mesh and write */
    PiercedStorage<double> data2(1, &patch_2D->getCells());
    data2.fill(0.0);
    {
        const PiercedStorage<mapping::Info> & mapper = mapobject.getInverseMapping();
        std::vector<double> vdata2(patch_2D->getInternalCount());
        count = 0;
        for (Cell & cell : patch_2D->getCells()){
            if (cell.isInterior()){
                long id = cell.getId();
                if (mapper[id].type == mapping::Type::TYPE_RENUMBERING){
                    data2[id] = data[mapper[id].mapped[0]];
                    vdata2[count] = data2[id];
                }
                else if (mapper[id].type == mapping::Type::TYPE_COARSENING){
                    data2[id] = 0.0;
                    int n = mapper[id].mapped.size();
                    for (long idd : mapper[id].mapped){
                        data2[id] += data[idd] / n;
                    }
                    vdata2[count] = data2[id];
                }
                else if (mapper[id].type == mapping::Type::TYPE_REFINEMENT){
                    data2[id] = data[mapper[id].mapped[0]];
                    vdata2[count] = data2[id];
                }
                count++;
            }
        }

        patch_2D->getVTK().setName("mesh_random.0");
        patch_2D->getVTK().addData("data", VTKFieldType::SCALAR, VTKLocation::CELL, vdata2);
        patch_2D->write();
        patch_2D->getVTK().setName("mesh_random.1");
        patch_2D->write();

    }

    /** Re-Map data on first mesh with inverse mapping and write */
    {
        const PiercedStorage<mapping::Info> & invmapper = mapobject.getMapping();
        PiercedStorage<double> data3(1, &patch_2D_original->getCells());
        std::vector<double> vdata3(patch_2D_original->getInternalCount());
        count = 0;
        for (Cell & cell : patch_2D_original->getCells()){
            if (cell.isInterior()){
                long id = cell.getId();
                if (invmapper[id].type == mapping::Type::TYPE_RENUMBERING){
                    data3[id] = data2[invmapper[id].mapped[0]];
                    vdata3[count] = data3[id];
                }
                else if (invmapper[id].type == mapping::Type::TYPE_COARSENING){
                    data3[id] = 0.0;
                    int n = invmapper[id].mapped.size();
                    for (long idd : invmapper[id].mapped){
                        data3[id] += data2[idd] / n;
                    }
                    vdata3[count] = data3[id];
                }
                else if (invmapper[id].type == mapping::Type::TYPE_REFINEMENT){
                    data3[id] = data2[invmapper[id].mapped[0]];
                    vdata3[count] = data3[id];
                }
                count++;
            }
        }

        patch_2D_original->getVTK().setName("mesh_original.1");
        patch_2D_original->getVTK().addData("data", VTKFieldType::SCALAR, VTKLocation::CELL, vdata3);
        patch_2D_original->write();

    }

}

/**
 * Main program.
 */

int main(int argc, char *argv[]) 
{
#if BITPIT_ENABLE_MPI==1
    MPI_Init(&argc,&argv);
#endif    

    /** Run the example **/
    try {
        run();
    } catch (const std::exception &exception) {
        log::cout() << exception.what();
        exit(1);
    }

#if BITPIT_ENABLE_MPI==1
    MPI_Finalize();
#endif

}
