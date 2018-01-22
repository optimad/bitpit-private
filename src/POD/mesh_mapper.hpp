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

#ifndef __BITPIT_MESH_MAPPER_HPP__
#define __BITPIT_MESH_MAPPER_HPP__

#if BITPIT_ENABLE_MPI
#    include <mpi.h>
#endif
#include <string>
#include <vector>
#include <unordered_map>

#include "bitpit_voloctree.hpp"

namespace bitpit {

namespace mapping
{

typedef adaption::Info Info;

enum Type {
    TYPE_UNKNOWN = 0,
    TYPE_REFINEMENT,
    TYPE_COARSENING,
    TYPE_RENUMBERING
};

enum Entity {
    ENTITY_UNKNOWN = -1,
    ENTITY_CELL,
};
struct mInfo
{
    mInfo()
    : type(TYPE_UNKNOWN), entity(ENTITY_UNKNOWN)
    {
    }

    mInfo(Type user_type, Entity user_entity)
    : type(user_type), entity(user_entity)
    {
    }

    Type type;
    Entity entity;
    std::vector<long> mapped;
# if BITPIT_ENABLE_MPI
    std::vector<int> rank;
# endif
};


}

class MeshMapper {

public:
# if BITPIT_ENABLE_MPI
    MeshMapper(MPI_Comm comm = MPI_COMM_WORLD);
# else
    MeshMapper();
# endif

    ~MeshMapper();
    MeshMapper(MeshMapper&& other) = default;

    void clear();
    void clearMapping();
    void clearInverseMapping();

    const bitpit::PiercedStorage<mapping::mInfo> & getMapping();
    const bitpit::PiercedStorage<mapping::mInfo> & getInverseMapping();

    void mapMeshes(bitpit::VolumeKernel * meshReference, bitpit::VolumeKernel * meshMapped, bool fillInv = false);

    void mappingAdaptionPreparare(const std::vector<adaption::Info> & infoAdapt, bool reference = true);
    void mappingAdaptionUpdate(const std::vector<adaption::Info> & infoAdapt, bool reference = true, bool fillInv = true);

protected:

#if BITPIT_ENABLE_MPI
    MPI_Comm                m_communicator; /**< MPI communicator */
    int                     m_rank;         /**< Local rank of process. */
    int                     m_nProcs;       /**< Number of processes. */
#endif

    VolumeKernel* m_referenceMesh;
    VolumeKernel* m_mappedMesh;

    PiercedStorage<mapping::mInfo> m_mapper;  /**< Mapping info for each cell of reference mesh.
                                                                  The mapping info is treated as a set of adaption info related to
                                                                  an adaption of the mapped mesh toward the reference mesh. */
    PiercedStorage<mapping::mInfo> m_invmapper;  /**< Inverse mapping info for each cell of mapped mesh. */


    std::unordered_map<long, mapping::mInfo> m_previousmapper;


    void _mapMeshes(bitpit::VolOctree * meshReference, bitpit::VolOctree * meshMapped, bool fillInv);
    void _mapMeshesSamePartition(bitpit::VolOctree * meshReference, bitpit::VolOctree * meshMapped, bool fillInv);

#if BITPIT_ENABLE_MPI
    void initializeCommunicator(MPI_Comm communicator);
    MPI_Comm getCommunicator() const;
    bool isCommunicatorSet() const;
    void freeCommunicator();
#endif

};

}
#endif
