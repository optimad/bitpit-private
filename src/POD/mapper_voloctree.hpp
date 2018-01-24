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

#ifndef __BITPIT_MAPPER_VOLOCTREE_HPP__
#define __BITPIT_MAPPER_VOLOCTREE_HPP__

#if BITPIT_ENABLE_MPI==1
#    include <mpi.h>
#endif
#include <string>
#include <vector>
#include <unordered_map>

#include "bitpit_voloctree.hpp"

namespace bitpit {

namespace mapping
{

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

/**
 * \class MappingInfo
 * \ingroup Ref
 *
 * \brief The MappingInfo is the structure to store info about mapping between elements.
 *
 */
struct MappingInfo
{

    /**
     * Default constructor.
     */
    MappingInfo()
    : type(TYPE_UNKNOWN), entity(ENTITY_UNKNOWN)
    {
    }

    /**
     * Custom constructor.
     * \param[in] user_type Type of mapping item
     * \param[in] user_entity Mapped entity
     */
    MappingInfo(Type user_type, Entity user_entity)
    : type(user_type), entity(user_entity)
    {
    }

    Type type;                  /**< Type of mapping item. */
    Entity entity;              /**< Mapped entity. */
    std::vector<long> mapped;   /**< ID of mapped elements. */
# if BITPIT_ENABLE_MPI==1
    std::vector<int> rank;      /**< Rank of mapped elements. */
# endif
};

}

class MapperVolOctree {

public:
# if BITPIT_ENABLE_MPI==1
    MapperVolOctree(VolOctree * meshReference, VolOctree * meshMapped, MPI_Comm comm = MPI_COMM_WORLD);
# else
    MapperVolOctree(VolOctree * meshReference, VolOctree * meshMapped);
# endif

    ~MapperVolOctree();

    /**
     * Copy constructor.
     */
    MapperVolOctree(MapperVolOctree&& other) = default;

    const bitpit::PiercedStorage<mapping::MappingInfo> & getMapping();
    const bitpit::PiercedStorage<mapping::MappingInfo> & getInverseMapping();

    void initialize(bool fillInv = false);

    void adaptionPrepare(const std::vector<adaption::Info> & infoAdapt, bool reference = true);
    void adaptionAlter(const std::vector<adaption::Info> & infoAdapt, bool reference = true, bool fillInv = false);
    void adaptionCleanup();

#if BITPIT_ENABLE_MPI==1
    void clearPartitionMapping();
    std::map<int, std::vector<long> > getReceivedMappedID();
    std::map<int, std::vector<long> > getSentMappedID();
#endif


protected:

    /**
     * \class OctantIR
     * \ingroup Ref
     *
     * \brief The OctantIR is the structure to store info about an octant in a mapper object.
     *
     */
    struct OctantIR{

        Octant octant;  /**< Octant object. */
        long ID;        /**< Octant local ID. */
        long globalID;  /**< Octant global ID. */
        int rank;       /**< Octant rank. */

        OctantIR(){};

        /**
         * Custom constructor.
         * \param[in] _octant Octant object
         * \param[in] _ID Local octant ID
         * \param[in] _globalID Global octant ID
         * \param[in] _rank Octant rank
         */
        OctantIR(Octant _octant, long _ID, long _globalID = -1, int _rank = -1)
        {
            octant = _octant;
            ID = _ID;
            globalID = _globalID;
            rank = _rank;
        };
    };

#if BITPIT_ENABLE_MPI==1

    /**
     * \class PartitionMapper
     * \ingroup Ref
     *
     * \brief The PartitionMapper is the structure to store info about the partitioning of two VolOctree mapped meshes.
     *
     */
    struct PartitionMapper {

        std::vector<OctantIR> list_rec_octantIR_before; /**< List of received Octants from lower rank processes. */
        std::vector<OctantIR> list_rec_octantIR_after;  /**< List of received Octants from higher rank processes. */

        std::unordered_map<int, std::unordered_map<long, OctantIR*> > map_rank_rec_octantIR;        /**< List of received Octants from other processes. */
        std::unordered_map<int, std::unordered_map<long, mapping::MappingInfo> > map_rank_invmapper;      /**< Inverse mapper terms related to processes different from the local rank. */
        std::unordered_map<int, std::unordered_map<long, mapping::MappingInfo> > map_rank_previousmapper; /**< Temporary inverse mapper for each processes used during a mesh adaptation. */

        std::vector<OctantIR> list_sent_octantIR;   /**< List of sent Octants to other processes. */

        std::vector<uint64_t> partitionFDReference;  /**< First descendant partitioning structure of reference mesh. */
        std::vector<uint64_t> partitionLDReference;  /**< Last descendant partitioning structure of reference mesh. */
        std::vector<uint64_t> partitionFDMapped;     /**< First descendant partitioning structure of mapped mesh. */
        std::vector<uint64_t> partitionLDMapped;     /**< Last descendant partitioning structure of mapped mesh. */


        /**
         * Default constructor.
         */
        PartitionMapper(){};

        /**
         * Clear members.
         */
        void clear(){
            list_rec_octantIR_before.clear();
            list_rec_octantIR_after.clear();
            list_sent_octantIR.clear();
            map_rank_rec_octantIR.clear();
            map_rank_invmapper.clear();
            map_rank_previousmapper.clear();
            partitionFDReference.clear();
            partitionLDReference.clear();
            partitionFDMapped.clear();
            partitionLDMapped.clear();
        }

        /**
         * Clear only list members.
         */
        void clearLists(){
            list_rec_octantIR_before.clear();
            list_rec_octantIR_after.clear();
            list_sent_octantIR.clear();
        }

    };
#endif

    VolOctree* m_referenceMesh;      /**< Pointer to reference mesh.*/
    VolOctree* m_mappedMesh;         /**< Pointer to mapped mesh.*/

    PiercedStorage<mapping::MappingInfo> m_mapper;  /**< Mapping info for each cell of reference mesh.
                                                                  The mapping info is treated as a set of adaption info related to
                                                                  an adaption of the mapped mesh toward the reference mesh. */
    PiercedStorage<mapping::MappingInfo> m_invmapper;  /**< Inverse mapping info for each cell of mapped mesh. */

    std::unordered_map<long, mapping::MappingInfo> m_previousmapper;   /**< Temporary mapping used during a mesh adaptation. */

#if BITPIT_ENABLE_MPI==1
    MPI_Comm                m_communicator; /**< MPI communicator */
    int                     m_rank;         /**< Local rank of process. */
    int                     m_nProcs;       /**< Number of processes. */
    PartitionMapper         m_partitionIR;  /**< Partitioning info structure. */

#endif

    void clear();
    void clearMapping();
    void clearInverseMapping();

    void _mapMeshes(bool fillInv);
    void _mapMeshesSamePartition(const std::vector<OctantIR> * octantsReference, const std::vector<OctantIR> * octantsMapped,
            bool fillInv, long & indRef);

    void _mappingAdaptionReferenceUpdate(const std::vector<adaption::Info> & infoAdapt, bool fillInv = false);
    void _mappingAdaptionMappedUpdate(const std::vector<adaption::Info> & infoAdapt);

#if BITPIT_ENABLE_MPI==1
    void initializeCommunicator(MPI_Comm communicator);
    MPI_Comm getCommunicator() const;
    bool isCommunicatorSet() const;
    void freeCommunicator();

    void clearPartitionMappingLists();
    bool _checkPartition(VolOctree * meshA, VolOctree * meshB);
    void _mapMeshPartitioned(bool fillInv);
    bool _recoverPartition();

    void _communicateInverseMapper(std::unordered_map<long, mapping::MappingInfo> & _invmapper, const std::vector<OctantIR> * octantsIRMapped);
    void _communicateInverseMapperBack();
    void _communicateMappedAdaptionInfo(const std::vector<adaption::Info> & infoAdapt, std::vector<adaption::Info> & MappingInfoAdapt);

#endif

};

}
#endif
