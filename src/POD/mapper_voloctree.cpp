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

#include "mapper_voloctree.hpp"

namespace bitpit {

/**
 * \class MapperVolOctree
 * \ingroup Ref
 *
 * \brief The MapperVolOctree is the class to map two meshes of class VolOctree.
 *
 * The MapperVolOctree allows to map meshes of class VolOctree.
 *
 */

/**
 * Creates a new MapperVolOctree object.
 *
 * \param[in] meshReference Pointer to reference mesh
 * \param[in] meshMapped Pointer to mapped mesh
 */
# if BITPIT_ENABLE_MPI==1
/**
 * \param[in] comm The MPI communicator used by the Ref object. MPI_COMM_WORLD is the default value.
 */
MapperVolOctree::MapperVolOctree(bitpit::VolOctree * meshReference, bitpit::VolOctree * meshMapped, MPI_Comm comm)
# else
MapperVolOctree::MapperVolOctree(bitpit::VolOctree * meshReference, bitpit::VolOctree * meshMapped)
# endif
: m_mapper(1)
{

    m_referenceMesh = meshReference;
    m_mappedMesh = meshMapped;

#if BITPIT_ENABLE_MPI==1
    m_communicator = MPI_COMM_NULL;
    initializeCommunicator(comm);
    MPI_Comm_size(m_communicator, &m_nProcs);
    MPI_Comm_rank(m_communicator, &m_rank);
#endif
}

/**
 * Destructor of MapperVolOctree
 */
MapperVolOctree::~MapperVolOctree()
{
}

/**
 * Clear mapping members
 */
void MapperVolOctree::clear()
{
    clearMapping();
    clearInverseMapping();
#if BITPIT_ENABLE_MPI==1
    clearPartitionMapping();
#endif
}

/**
 * Clear direct mapping
 */
void MapperVolOctree::clearMapping()
{
    if (m_mapper.getKernel() != nullptr)
        m_mapper.unsetKernel(true);
}

/**
 * Clear inverse mapping
 */
void MapperVolOctree::clearInverseMapping()
{
    if (m_invmapper.getKernel() != nullptr)
        m_invmapper.unsetKernel(false);
}

#if BITPIT_ENABLE_MPI==1
/**
 * Clear partition mapping members
 */
void MapperVolOctree::clearPartitionMapping()
{
    m_partitionIR.clear();
}

/**
 * Clear partition mapping lists
 */
void MapperVolOctree::clearPartitionMappingLists()
{
    m_partitionIR.clearLists();
}
#endif

/**
 * Get direct mapping
 */
const bitpit::PiercedStorage<mapping::MappingInfo> & MapperVolOctree::getMapping()
{
    return m_mapper;
}

/**
 * Get inverse mapping
 */
const bitpit::PiercedStorage<mapping::MappingInfo> & MapperVolOctree::getInverseMapping()
{
    return m_invmapper;
}

/**
 * Initialize the mapper of the mapped mesh on the reference mesh.
 * The two VolOctree meshes are set in the constructor of the object.
 *
 * \param[in] fillInv True if inverse mapper has to be computed
 */
void MapperVolOctree::initialize(bool fillInv)
{
    clear();
    _mapMeshes(fillInv);
}

/**
 * Prepare the mapper for an adaption of the reference OR the mapped mesh.
 * The adaptation of only one mesh is allowed.
 *
 * \param[in] infoAdapt Structure of adaptation info given by adaptionPrepare method of a patch
 * \param[in] reference True if the reference mesh will be adapted, false if the mapped one
 */
void MapperVolOctree::adaptionPrepare(const std::vector<adaption::Info> & infoAdapt, bool reference)
{
    m_previousmapper.clear();
    PiercedStorage<mapping::MappingInfo>* pmapper;

    if (reference)
        pmapper = &m_mapper;
    else
        pmapper = &m_invmapper;

    m_previousmapper.clear();
#if BITPIT_ENABLE_MPI==1
    m_partitionIR.map_rank_previousmapper.clear();
#endif
    for (const adaption::Info & info : infoAdapt){
        if (info.type == adaption::Type::TYPE_PARTITION_SEND ||
                info.type == adaption::Type::TYPE_PARTITION_RECV ||
                info.type == adaption::Type::TYPE_PARTITION_NOTICE){
            throw std::runtime_error ("Mapper: type of adation not supported : " + std::to_string(info.type));
        }
        else{
            if (info.type != adaption::Type::TYPE_DELETION &&
                    info.type != adaption::Type::TYPE_CREATION){
                for (const long & id : info.previous){
                    m_previousmapper[id] = (*pmapper)[id];
                }
            }
        }
    }
}

/**
 * Update the mapper after an adaption of the reference OR the mapped mesh.
 *
 * \param[in] infoAdapt Structure of adaptation info given by adaptionAlter method of a patch
 * \param[in] reference True if the reference mesh is adapted, false if the mapped one
 * \param[in] fillInv True if the inverse mapped was filled during the mesh mapper computing
 * (Note: if the adapted mesh is the mapped one, the inverse mapper has to be filled inevitably).
 */
void MapperVolOctree::adaptionAlter(const std::vector<adaption::Info> & infoAdapt, bool reference, bool fillInv)
{
    if (reference)
        _mappingAdaptionReferenceUpdate(infoAdapt, fillInv);
    else
        _mappingAdaptionMappedUpdate(infoAdapt);
}

/**
 * Clear the mapper internal structures.
 */
void MapperVolOctree::adaptionCleanup(){
    m_previousmapper.clear();
}

/**
 * Update the mapper after an adaption of the reference mesh.
 *
 * \param[in] infoAdapt Structure of adaptation info given by adaptionAlter method of a patch
 * \param[in] fillInv True if the inverse mapped was filled during the mesh mapper computing
 */
void MapperVolOctree::_mappingAdaptionReferenceUpdate(const std::vector<adaption::Info> & infoAdapt, bool fillInv)
{

    VolOctree* meshAdapted;
    PiercedStorage<mapping::MappingInfo>* mapperAdapted;
    VolOctree* meshMapped;
    PiercedStorage<mapping::MappingInfo>* mapperMapped;

    meshAdapted = m_referenceMesh;
    mapperAdapted = &m_mapper;
    meshMapped = m_mappedMesh;
    mapperMapped = &m_invmapper;

    bool checkPart = true;
    bool changedPartition = false;
#if BITPIT_ENABLE_MPI==1
    checkPart = checkPartition();
    if (!checkPart)
        changedPartition = _recoverPartition();
#endif

    if (!changedPartition){

        for (const adaption::Info & info : infoAdapt){

            if (info.type == adaption::Type::TYPE_PARTITION_SEND ||
                    info.type == adaption::Type::TYPE_PARTITION_RECV ||
                    info.type == adaption::Type::TYPE_PARTITION_NOTICE){
                throw std::runtime_error ("Mapper: type of adation not supported : " + std::to_string(info.type));
            }
            else if (info.type == adaption::Type::TYPE_DELETION ||
                    info.type == adaption::Type::TYPE_CREATION){
                //Do nothing
            }
            else{

                if (fillInv){
                    //Erase previous id in mapperMapped elements
                    for (long idprevious : info.previous){
                        std::vector<long> _mapped = m_previousmapper[idprevious].mapped;
                        std::size_t imapped = 0;
                        for (long idp : _mapped){
#if BITPIT_ENABLE_MPI==1
                            if (checkPart
                                    || m_previousmapper[idprevious].rank[imapped] == m_rank
                            ){
#endif
                                std::vector<long>::iterator it = std::find((*mapperMapped)[idp].mapped.begin(), (*mapperMapped)[idp].mapped.end(), idprevious);
                                int dist = std::distance((*mapperMapped)[idp].mapped.begin(), it);
                                (*mapperMapped)[idp].mapped.erase(it);
#if BITPIT_ENABLE_MPI==1
                                (*mapperMapped)[idp].rank.erase((*mapperMapped)[idp].rank.begin()+dist);
                            }
                            else{
                                mapping::MappingInfo & MappingInfo = m_partitionIR.map_rank_invmapper[m_previousmapper[idprevious].rank[imapped]][idp];
                                std::vector<long>::iterator it = std::find(MappingInfo.mapped.begin(), MappingInfo.mapped.end(), idprevious);
                                int dist = std::distance(MappingInfo.mapped.begin(), it);
                                MappingInfo.mapped.erase(it);
                                MappingInfo.rank.erase(MappingInfo.rank.begin()+dist);
                            }
#endif
                            imapped++;
                        }
                    }
                }

                switch(info.type){

                //RENUMBERING
                case adaption::Type::TYPE_RENUMBERING:
                {
                    long id = info.current[0];
                    std::vector<long> _mapped;
                    (*mapperAdapted)[id].mapped.clear();
                    (*mapperAdapted)[id].type = m_previousmapper[info.previous[0]].type;
                    (*mapperAdapted)[id].entity = mapping::Entity::ENTITY_CELL;
                    (*mapperAdapted)[id].mapped = m_previousmapper[info.previous[0]].mapped;
#if BITPIT_ENABLE_MPI==1
                    (*mapperAdapted)[id].rank = m_previousmapper[info.previous[0]].rank;
#endif
                    _mapped = (*mapperAdapted)[id].mapped;

                    if (fillInv){
                        std::size_t imapped = 0;
                        for (long idp : (*mapperAdapted)[id].mapped){
#if BITPIT_ENABLE_MPI==1
                            if (checkPart
                                    || (*mapperAdapted)[id].rank[imapped] == m_rank
                            ){
#endif
                                (*mapperMapped)[idp].mapped.push_back(id);
#if BITPIT_ENABLE_MPI==1
                                (*mapperMapped)[idp].rank.push_back(info.rank);
                            }
                            else{
                                mapping::MappingInfo & MappingInfo = m_partitionIR.map_rank_invmapper[(*mapperAdapted)[id].rank[imapped]][idp];
                                MappingInfo.mapped.push_back(id);
                                MappingInfo.rank.push_back(info.rank);
                            }
#endif
                            imapped++;
                        }
                    }
                }
                break;

                //REFINEMENT & COARSENING
                case adaption::Type::TYPE_REFINEMENT:
                case adaption::Type::TYPE_COARSENING:

                    //Clear current mapper
                    for (long id : info.current){
                        (*mapperAdapted)[id].mapped.clear();
#if BITPIT_ENABLE_MPI==1
                        (*mapperAdapted)[id].rank.clear();
#endif

                        (*mapperAdapted)[id].entity = mapping::Entity::ENTITY_CELL;
                    }

                    for (long idprevious : info.previous){
                        std::vector<long> _mapped;
                        std::vector<int> _rankmapped;
                        _mapped =  m_previousmapper[idprevious].mapped;
#if BITPIT_ENABLE_MPI==1
                        _rankmapped =  m_previousmapper[idprevious].rank;
#endif
                        std::size_t imapped = 0;
                        for (long idpmapped : _mapped){
                            VolOctree::OctantInfo oinfoprev;
                            Octant * octm ;
                            uint64_t mortonmapped;
                            uint64_t mortonlastdescmapped;
                            uint8_t levelmapped;
                            int rankmapped;

#if BITPIT_ENABLE_MPI==1
                            if (checkPart
                                    || _rankmapped[imapped] == m_rank
                            ){
#endif
                                oinfoprev = meshMapped->getCellOctant(idpmapped);
                                octm = meshMapped->getOctantPointer(oinfoprev);
                                mortonmapped = meshMapped->getTree().getMorton(octm);
                                mortonlastdescmapped = meshMapped->getTree().getLastDescMorton(octm);
                                levelmapped = meshMapped->getCellLevel(idpmapped);
#if BITPIT_ENABLE_MPI==1
                                rankmapped = _rankmapped[imapped];
                            }
                            else{
                                OctantIR* poct = m_partitionIR.map_rank_rec_octantIR[_rankmapped[imapped]][idpmapped];
                                octm = &poct->octant;
                                mortonmapped = meshMapped->getTree().getMorton(&poct->octant);
                                mortonlastdescmapped = meshMapped->getTree().getLastDescMorton(&poct->octant);
                                levelmapped = meshMapped->getTree().getLevel(&poct->octant);
                                rankmapped = _rankmapped[imapped];
                            }
#endif

                            for (long id : info.current){
                                //Retrieve level of current cell
                                uint8_t level;
                                uint64_t morton, mortonlastdesc;
                                level = meshAdapted->getCellLevel(id);
                                VolOctree::OctantInfo octinfo = meshAdapted->getCellOctant(id);
                                Octant * oct = meshAdapted->getOctantPointer(octinfo);
                                morton = meshAdapted->getTree().getMorton(oct);
                                mortonlastdesc = meshAdapted->getTree().getLastDescMorton(oct);

                                bool checkmorton = false;

                                checkmorton |= (morton >= mortonmapped && mortonlastdesc <= mortonlastdescmapped);
                                checkmorton |= (morton <= mortonmapped && mortonlastdesc >= mortonlastdescmapped);

                                if (checkmorton){

                                    if (level == levelmapped){
                                        (*mapperAdapted)[id].type = mapping::Type::TYPE_RENUMBERING;
                                        (*mapperAdapted)[id].mapped.push_back(idpmapped);
#if BITPIT_ENABLE_MPI==1
                                        (*mapperAdapted)[id].rank.push_back(rankmapped);

                                        if (fillInv){
                                            if (checkPart
                                                    || rankmapped == m_rank
                                            ){
#endif
                                                (*mapperMapped)[idpmapped].type = mapping::Type::TYPE_RENUMBERING;
                                                (*mapperMapped)[idpmapped].mapped.clear();
                                                (*mapperMapped)[idpmapped].mapped.push_back(id);
#if BITPIT_ENABLE_MPI==1
                                                (*mapperMapped)[idpmapped].rank.push_back(info.rank);
                                            }
                                            else{
                                                m_partitionIR.map_rank_invmapper[rankmapped][idpmapped].type = mapping::Type::TYPE_RENUMBERING;
                                                m_partitionIR.map_rank_invmapper[rankmapped][idpmapped].mapped.push_back(id);
                                                m_partitionIR.map_rank_invmapper[rankmapped][idpmapped].rank.push_back(info.rank);
                                            }
                                        }
#endif
                                    }
                                    else if (level > levelmapped){
                                        (*mapperAdapted)[id].type = mapping::Type::TYPE_REFINEMENT;
                                        (*mapperAdapted)[id].mapped.push_back(idpmapped);
#if BITPIT_ENABLE_MPI==1
                                        (*mapperAdapted)[id].rank.push_back(rankmapped);

                                        if (fillInv){
                                            if (checkPart
                                                    || rankmapped == m_rank
                                            ){
#endif
                                                (*mapperMapped)[idpmapped].mapped.push_back(id);
                                                (*mapperMapped)[idpmapped].type = mapping::Type::TYPE_COARSENING;
#if BITPIT_ENABLE_MPI==1
                                                (*mapperMapped)[idpmapped].rank.push_back(info.rank);
                                            }
                                            else{
                                                m_partitionIR.map_rank_invmapper[rankmapped][idpmapped].type = mapping::Type::TYPE_COARSENING;
                                                m_partitionIR.map_rank_invmapper[rankmapped][idpmapped].mapped.push_back(id);
                                                m_partitionIR.map_rank_invmapper[rankmapped][idpmapped].rank.push_back(info.rank);
                                            }
                                        }
#endif
                                    }
                                    else if (level < levelmapped){
                                        (*mapperAdapted)[id].type = mapping::Type::TYPE_COARSENING;
                                        (*mapperAdapted)[id].mapped.push_back(idpmapped);
#if BITPIT_ENABLE_MPI==1
                                        (*mapperAdapted)[id].rank.push_back(rankmapped);

                                        if (fillInv){
                                            if (checkPart
                                                    || rankmapped == m_rank
                                            ){
#endif
                                                (*mapperMapped)[idpmapped].type = mapping::Type::TYPE_REFINEMENT;
                                                (*mapperMapped)[idpmapped].mapped.push_back(id);
#if BITPIT_ENABLE_MPI==1
                                                (*mapperMapped)[idpmapped].rank.push_back(info.rank);
                                            }
                                            else{
                                                m_partitionIR.map_rank_invmapper[rankmapped][idpmapped].type = mapping::Type::TYPE_REFINEMENT;
                                                m_partitionIR.map_rank_invmapper[rankmapped][idpmapped].mapped.push_back(id);
                                                m_partitionIR.map_rank_invmapper[rankmapped][idpmapped].rank.push_back(info.rank);
                                            }
                                        }
#endif

                                    } //end level < levelreference
                                } //end checkmorton
                            } //end id current

                            imapped++;

                        }//idmapped
                    } //idprevious
                    break;

                }//end switch
            }//if deletion
        }//end info


#if BITPIT_ENABLE_MPI==1
        if (fillInv)
            _communicateInverseMapperBack();
#endif

    }
    else{

        initialize(fillInv);

    }//endif changedPartition
}

/**
 * Update the mapper after an adaption of the reference OR the mapped mesh.
 *
 * \param[in] infoAdapt Structure of adaptation info given by adaptionAlter method of a patch
 */
void MapperVolOctree::_mappingAdaptionMappedUpdate(const std::vector<adaption::Info> & infoAdapt)
{

    VolOctree* meshAdapted;
    PiercedStorage<mapping::MappingInfo>* mapperAdapted;
    VolOctree* meshReference;
    PiercedStorage<mapping::MappingInfo>* mapperReference;

    meshAdapted = m_mappedMesh;
    mapperAdapted = &m_invmapper;
    meshReference = m_referenceMesh;
    mapperReference = &m_mapper;

    bool checkPart = true;
    bool changedPartition = false;
#if BITPIT_ENABLE_MPI==1
    checkPart = checkPartition();
    if (!checkPart)
        changedPartition = _recoverPartition();
#endif

    if (!m_invmapper.getKernel())
        throw std::runtime_error("update mapper with adaption of mapped mesh possible only if inverse mapper is filled");

    if (!changedPartition){

        std::vector<adaption::Info> MappingInfoAdapt;
#if BITPIT_ENABLE_MPI==1
        if (!checkPart)
            _communicateMappedAdaptionInfo(infoAdapt, MappingInfoAdapt);
        else
#endif
            MappingInfoAdapt = infoAdapt;

        for (const adaption::Info & info : MappingInfoAdapt){

            if (info.type == adaption::Type::TYPE_PARTITION_SEND ||
                    info.type == adaption::Type::TYPE_PARTITION_RECV ||
                    info.type == adaption::Type::TYPE_PARTITION_NOTICE){
                throw std::runtime_error ("Mapper: type of adation not supported : " + std::to_string(info.type));
            }
            else if (info.type == adaption::Type::TYPE_DELETION ||
                    info.type == adaption::Type::TYPE_CREATION){
                //Do nothing
            }
            else{

                //Erase previous id in mapperReference elements
                for (long idprevious : info.previous){
                    std::vector<long> _mapped;
#if BITPIT_ENABLE_MPI==1
                    if (checkPart
                            || info.rank == m_rank
                    ){
#endif
                        _mapped =  m_previousmapper[idprevious].mapped;
#if BITPIT_ENABLE_MPI==1
                    }
                    else{
                        _mapped = m_partitionIR.map_rank_previousmapper[info.rank][idprevious].mapped;
                    }
#endif
                    for (long idp : _mapped){
                        std::vector<long>::iterator it = std::find((*mapperReference)[idp].mapped.begin(), (*mapperReference)[idp].mapped.end(), idprevious);
                        int dist = std::distance((*mapperReference)[idp].mapped.begin(), it);
                        (*mapperReference)[idp].mapped.erase(it);
#if BITPIT_ENABLE_MPI==1
                        (*mapperReference)[idp].rank.erase((*mapperReference)[idp].rank.begin()+dist);
#endif
                    }
                }


                //RENUMBERING
                if (info.type == adaption::Type::TYPE_RENUMBERING){
                    long id = info.current[0];
                    std::vector<long> _mapped;
#if BITPIT_ENABLE_MPI==1
                    if (checkPart
                            || info.rank == m_rank
                    ){
#endif
                        (*mapperAdapted)[id].mapped.clear();
                        (*mapperAdapted)[id].type = m_previousmapper[info.previous[0]].type;
                        (*mapperAdapted)[id].entity = mapping::Entity::ENTITY_CELL;
                        (*mapperAdapted)[id].mapped = m_previousmapper[info.previous[0]].mapped;
                        _mapped = (*mapperAdapted)[id].mapped;
#if BITPIT_ENABLE_MPI==1
                        (*mapperAdapted)[id].rank = m_previousmapper[info.previous[0]].rank;
                    }
                    else{
                        for (const long & _idp : info.previous){
                            m_partitionIR.map_rank_invmapper[info.rank].erase(_idp);
                        }

                        m_partitionIR.map_rank_invmapper[info.rank][id].mapped.clear();
                        m_partitionIR.map_rank_invmapper[info.rank][id].type = m_partitionIR.map_rank_previousmapper[info.rank][info.previous[0]].type;
                        m_partitionIR.map_rank_invmapper[info.rank][id].entity = mapping::Entity::ENTITY_CELL;
                        m_partitionIR.map_rank_invmapper[info.rank][id].mapped = m_partitionIR.map_rank_previousmapper[info.rank][info.previous[0]].mapped;
                        m_partitionIR.map_rank_invmapper[info.rank][id].rank = m_partitionIR.map_rank_previousmapper[info.rank][info.previous[0]].rank;
                        _mapped = m_partitionIR.map_rank_invmapper[info.rank][id].mapped;
                    }
#endif
                    std::size_t imapped = 0;
                    for (long idp : _mapped){
                        (*mapperReference)[idp].mapped.push_back(id);
#if BITPIT_ENABLE_MPI==1
                        (*mapperReference)[idp].rank.push_back(info.rank);
#endif
                    }
                }

                //REFINEMENT & COARSENING
                if (info.type == adaption::Type::TYPE_REFINEMENT || info.type == adaption::Type::TYPE_COARSENING){

                    //Clear current mapper
                    for (long id : info.current){

#if BITPIT_ENABLE_MPI==1
                        if (checkPart
                                || info.rank == m_rank
                        ){
#endif
                            (*mapperAdapted)[id].mapped.clear();
                            (*mapperAdapted)[id].entity = mapping::Entity::ENTITY_CELL;
#if BITPIT_ENABLE_MPI==1
                            (*mapperAdapted)[id].rank.clear();
                        }
                        else{
                            m_partitionIR.map_rank_invmapper[info.rank][id].mapped.clear();
                            m_partitionIR.map_rank_invmapper[info.rank][id].rank.clear();
                            m_partitionIR.map_rank_invmapper[info.rank][id].entity = mapping::Entity::ENTITY_CELL;
                        }
#endif
                    }

                    for (long idprevious : info.previous){
                        std::vector<long> _mapped;
                        std::vector<int> _rankmapped;
#if BITPIT_ENABLE_MPI==1
                        if (checkPart
                                || info.rank == m_rank
                        ){
#endif
                            _mapped =  m_previousmapper[idprevious].mapped;
#if BITPIT_ENABLE_MPI==1
                            _rankmapped =  m_previousmapper[idprevious].rank;
                        }
                        else{
                            _mapped = m_partitionIR.map_rank_previousmapper[info.rank][idprevious].mapped;
                            _rankmapped = m_partitionIR.map_rank_previousmapper[info.rank][idprevious].rank;
                        }
#endif
                        std::size_t icount = 0;
                        for (long idpmapped : _mapped){

                            VolOctree::OctantInfo oinfoprev = meshReference->getCellOctant(idpmapped);
                            Octant * oct = meshReference->getOctantPointer(oinfoprev);
                            uint64_t mortonreference = meshReference->getTree().getMorton(oct);
                            uint64_t mortonlastdescreference = meshReference->getTree().getLastDescMorton(oct);
                            uint8_t levelreference = meshReference->getCellLevel(idpmapped);
#if BITPIT_ENABLE_MPI==1
                            //IT'S ALWAYS = m_rank?!...
                            int rankmapped = _rankmapped[icount];
#endif
                            for (long id : info.current){

                                //Retrieve level of current cell
                                uint8_t level;
                                uint64_t morton, mortonlastdesc;
#if BITPIT_ENABLE_MPI==1
                                if (checkPart
                                        || info.rank == m_rank
                                ){
#endif
                                    level = meshAdapted->getCellLevel(id);
                                    VolOctree::OctantInfo octinfo = meshAdapted->getCellOctant(id);
                                    Octant * oct = meshAdapted->getOctantPointer(octinfo);
                                    morton = meshAdapted->getTree().getMorton(oct);
                                    mortonlastdesc = meshAdapted->getTree().getLastDescMorton(oct);
#if BITPIT_ENABLE_MPI==1
                                }
                                else{
                                    OctantIR* poct = m_partitionIR.map_rank_rec_octantIR[info.rank][id];
                                    level = poct->octant.getLevel();
                                    morton = meshAdapted->getTree().getMorton(&poct->octant);
                                    mortonlastdesc = meshAdapted->getTree().getLastDescMorton(&poct->octant);
                                }
#endif

                                bool checkmorton = false;

                                checkmorton |= (morton >= mortonreference && mortonlastdesc <= mortonlastdescreference);
                                checkmorton |= (morton <= mortonreference && mortonlastdesc >= mortonlastdescreference);

                                if (checkmorton){

                                    if (level == levelreference){
#if BITPIT_ENABLE_MPI==1
                                        if (checkPart
                                                || info.rank == m_rank
                                        ){
#endif
                                            (*mapperAdapted)[id].type = mapping::Type::TYPE_RENUMBERING;
                                            (*mapperAdapted)[id].mapped.push_back(idpmapped);
#if BITPIT_ENABLE_MPI==1
                                            (*mapperAdapted)[id].rank.push_back(info.rank);
                                        }
                                        else{
                                            m_partitionIR.map_rank_invmapper[info.rank][id].type = mapping::Type::TYPE_RENUMBERING;
                                            m_partitionIR.map_rank_invmapper[info.rank][id].mapped.push_back(idpmapped);
                                            m_partitionIR.map_rank_invmapper[info.rank][id].rank.push_back(rankmapped);
                                        }
                                        (*mapperReference)[idpmapped].rank.push_back(info.rank);
#endif
                                        (*mapperReference)[idpmapped].type = mapping::Type::TYPE_RENUMBERING;
                                        (*mapperReference)[idpmapped].mapped.clear();
                                        (*mapperReference)[idpmapped].mapped.push_back(id);
                                    }
                                    else if (level > levelreference){
#if BITPIT_ENABLE_MPI==1
                                        if (checkPart
                                                || info.rank == m_rank
                                        ){
#endif
                                            (*mapperAdapted)[id].type = mapping::Type::TYPE_REFINEMENT;
                                            (*mapperAdapted)[id].mapped.push_back(idpmapped);
#if BITPIT_ENABLE_MPI==1
                                            (*mapperAdapted)[id].rank.push_back(rankmapped);
                                        }
                                        else{
                                            m_partitionIR.map_rank_invmapper[info.rank][id].type = mapping::Type::TYPE_REFINEMENT;
                                            m_partitionIR.map_rank_invmapper[info.rank][id].mapped.push_back(idpmapped);
                                            m_partitionIR.map_rank_invmapper[info.rank][id].rank.push_back(rankmapped);
                                        }
                                        (*mapperReference)[idpmapped].rank.push_back(info.rank);
#endif
                                        (*mapperReference)[idpmapped].mapped.push_back(id);
                                        (*mapperReference)[idpmapped].type = mapping::Type::TYPE_COARSENING;
                                    }
                                    else if (level < levelreference){

#if BITPIT_ENABLE_MPI==1
                                        if (checkPart
                                                || info.rank == m_rank
                                        ){
#endif
                                            (*mapperAdapted)[id].type = mapping::Type::TYPE_COARSENING;
                                            (*mapperAdapted)[id].mapped.push_back(idpmapped);
#if BITPIT_ENABLE_MPI==1
                                            (*mapperAdapted)[id].rank.push_back(rankmapped);
                                        }
                                        else{
                                            m_partitionIR.map_rank_invmapper[info.rank][id].type = mapping::Type::TYPE_COARSENING;
                                            m_partitionIR.map_rank_invmapper[info.rank][id].mapped.push_back(idpmapped);
                                            m_partitionIR.map_rank_invmapper[info.rank][id].rank.push_back(rankmapped);
                                        }
                                        (*mapperReference)[idpmapped].rank.push_back(info.rank);
#endif
                                        (*mapperReference)[idpmapped].type = mapping::Type::TYPE_REFINEMENT;
                                        (*mapperReference)[idpmapped].mapped.push_back(id);

                                    } //end level < levelreference
                                } //end checkmorton
                            } //end id current
                            icount++;
                        }//idmapped
                    } //idprevious
                }
            }//if deletion
        }//end info

#if BITPIT_ENABLE_MPI==1
        _communicateInverseMapperBack();
#endif

    }
    else{

        initialize(true);

    }//endif changedPartition
}


#if BITPIT_ENABLE_MPI==1

/**
 * Get the list of the octants of the partitions, different from the local one, of the mapped mesh overlapped to the local partition of the reference mesh.
 *
 * return Map with for each rank (key) the list of the octants (argument) of the mapped mesh overlapped with the local partition of the reference mesh.
 */
std::map<int, std::vector<long> > MapperVolOctree::getReceivedMappedID()
{
    std::map<int, std::vector<long> > received;

    for (OctantIR octir : m_partitionIR.list_rec_octantIR_before){
        received[octir.rank].push_back(octir.ID);
    }
    for (OctantIR octir : m_partitionIR.list_rec_octantIR_after){
        received[octir.rank].push_back(octir.ID);
    }

    return received;

}

/**
 * Get the list of the octants of the local partition of the reference mesh overlapped
 * to a different partition of the mapped mesh.
 *
 * return Map with for each rank (key) the list of the octants (argument) of the reference mesh overlapped with the partition (rank) of the mapped mesh.
 */
std::map<int, std::vector<long> > MapperVolOctree::getSentReferenceID()
{
    std::map<int, std::vector<long> > sent;

    {
        // recover id to be rec/sent
        std::map<int, std::set<long> > rankIDsend;

        for (Cell & cell : m_referenceMesh->getCells()){
            long ID = cell.getId();
            auto info = m_mapper[ID];
            int i = 0;
            for (int rank : info.rank){
                if (rank != m_referenceMesh->getRank()){
                    rankIDsend[rank].insert(ID);
                }
                i++;
            }
        }

        for (int irank=0; irank<m_nProcs; irank++){
            sent[irank].reserve(rankIDsend[irank].size());
        }

        for (std::map<int, std::set<long> >::iterator it=rankIDsend.begin(); it!=rankIDsend.end(); it++)
        {
            for (long id : it->second)
                sent[it->first].push_back(id);
        }
    }

    return sent;

}

/**
 * Get the list of the octants of the local partition of the mapped mesh overlapped to a different partition of the reference mesh.
 *
 * return Map with for each rank (key) the list of the octants (argument) of the mapped mesh overlapped with the partition (rank) of the reference mesh.
 */
std::map<int, std::vector<long> > MapperVolOctree::getSentMappedID()
{
    std::map<int, std::vector<long> > sent;

    for (OctantIR octir : m_partitionIR.list_sent_octantIR){
        sent[octir.rank].push_back(octir.ID);
    }

    return sent;

}

#endif


/**
 * Map an input VolOctree mesh on a VolOctree reference mesh already set in constructor.
 * Requirement : the meshes have to be defined on the same identical domain.
 *
 * \param[in] fillInv If true even the inverse mapping (reference mesh to input mesh) is filled.
 */
void MapperVolOctree::_mapMeshes(bool fillInv)
{

    if ( (m_referenceMesh->getLength() != m_mappedMesh->getLength()) || (m_referenceMesh->getOrigin() != m_mappedMesh->getOrigin()) )
        throw std::runtime_error ("mesh mapper: different domain of VolOctree meshes not allowed.");

#if BITPIT_ENABLE_MPI==1
    clearPartitionMapping();
#endif

    m_mapper.setDynamicKernel(&m_referenceMesh->getCells(), PiercedSyncMaster::SyncMode::SYNC_MODE_JOURNALED);
    if (fillInv)
        m_invmapper.setDynamicKernel(&m_mappedMesh->getCells(), PiercedSyncMaster::SyncMode::SYNC_MODE_JOURNALED);
    else
        clearInverseMapping();


    bool checkPart = true;

#if BITPIT_ENABLE_MPI==1

    if (!(m_referenceMesh->isPartitioned())){
#endif

        //        _mapMeshesSamePartition(m_referenceMesh, m_mappedMesh, fillInv);
        long indRef = 0;
        _mapMeshesSamePartition(nullptr, nullptr, fillInv, indRef);

#if BITPIT_ENABLE_MPI==1
    }
    else{

        checkPart = checkPartition();

        if (checkPart){
            long indRef = 0;
            _mapMeshesSamePartition(nullptr, nullptr, fillInv, indRef);
        }
        else{
            _mapMeshPartitioned(fillInv);
        }

    }
#endif

}

/**
 * Map an input list of octants of a VolOctree mesh on a list of octants of a VolOctree reference mesh.
 * Requirements : the meshes have to be defined on the same identical domain;
 *                the lists are considered as identically parallel partitioned,
 *                i.e. the first and the last descendant octants of the input list are
 *                contained in the covered domain of the reference list.
 *
 * \param[in] octantsIRReference Pointer to reference octant list (if =nullptr the whole local mesh is used)
 * \param[in] octantsIRMapped Pointer to reference mesh (if =nullptr the whole local mesh is used)
 * \param[in] fillInv If true even the inverse mapping (reference mesh to input mesh) is filled [NOTE: forced to false in this version!].
 * \param[in, out] indRef Starting octant (ending octant in output) of reference mesh to compute the mapping (default = 0).
 */
void MapperVolOctree::_mapMeshesSamePartition(const std::vector<OctantIR> * octantsIRReference, const std::vector<OctantIR> * octantsIRMapped,
        bool fillInv, long & indRef)
{

    //VolOctree specialization
    bitpit::VolOctree* meshReference = m_referenceMesh;
    bitpit::VolOctree* meshMapped = m_mappedMesh;

    //Fill IR with meshes if list pointer is null
    std::vector<OctantIR> tempOctantsIRReference;
    std::vector<OctantIR> tempOctantsIRMapped;
    bool localMapped = false;
    if (octantsIRReference == nullptr){
        long n = meshReference->getInternalCount();
        const std::vector<Octant>* octants = meshReference->getTree().getInternalOctants();
        for (long i = 0; i < n; i++){
            VolOctree::OctantInfo octinfoRef(i, true);
            long idRef = meshReference->getOctantId(octinfoRef);
            OctantIR val(octants->at(i), idRef, idRef);
#if BITPIT_ENABLE_MPI==1
            val.rank = m_rank;
#endif
            tempOctantsIRReference.push_back(val);
        }
        octantsIRReference = &tempOctantsIRReference;
    }
    if (octantsIRMapped == nullptr){
        long n = meshMapped->getInternalCount();
        const std::vector<Octant>* octants = meshMapped->getTree().getInternalOctants();
        for (long i = 0; i < n; i++){
            VolOctree::OctantInfo octinfoMap(i, true);
            long idMap = meshMapped->getOctantId(octinfoMap);
            OctantIR val(octants->at(i), idMap, idMap);
#if BITPIT_ENABLE_MPI==1
            val.rank = m_rank;
#endif
            tempOctantsIRMapped.push_back(val);
        }
        octantsIRMapped = &tempOctantsIRMapped;
        localMapped = true;
    }

    //Define a map for inverse mapper structure.
    //If serial the map is transfer to inverse mapper directly.
    //If mapped mesh is parallel and partitioned the map object
    //is communicated at the end of the procedure.
    std::unordered_map<long, mapping::MappingInfo> _invmapper;

    long nRef  = octantsIRReference->size();
    long nMap  = octantsIRMapped->size();

    long indMap = 0;

    if (indRef < nRef && indMap < nMap)
    {
        //Place indMap at first last desc morton greater than first morton of reference octants
        uint64_t mortonMaplastdesc = meshMapped->getTree().getLastDescMorton(&octantsIRMapped->at(indMap).octant);
        uint64_t morton = meshReference->getTree().getMorton(&octantsIRReference->at(indRef).octant);
        while (mortonMaplastdesc < morton && indMap < nMap){
            indMap++;
            mortonMaplastdesc = meshMapped->getTree().getLastDescMorton(&octantsIRMapped->at(indMap).octant);
        }
    }

    while (indRef < nRef && indMap < nMap){
        long idRef = octantsIRReference->at(indRef).ID;
        long idMap = octantsIRMapped->at(indMap).ID;
        long gidMap = octantsIRMapped->at(indMap).globalID;
        int rank = octantsIRMapped->at(indMap).rank;

        m_mapper[idRef].entity = bitpit::mapping::Entity::ENTITY_CELL;
        if (fillInv){
            _invmapper[gidMap].entity = bitpit::mapping::Entity::ENTITY_CELL;
        }

        if (octantsIRMapped->at(indMap).octant.getLevel() == octantsIRReference->at(indRef).octant.getLevel()){
            m_mapper[idRef].mapped.push_back(idMap);
            m_mapper[idRef].type = bitpit::mapping::Type::TYPE_RENUMBERING;
#if BITPIT_ENABLE_MPI==1
            m_mapper[idRef].rank.push_back(rank);
#endif
            if (fillInv){
                _invmapper[gidMap].mapped.push_back(idRef);
                _invmapper[gidMap].type = bitpit::mapping::Type::TYPE_RENUMBERING;
#if BITPIT_ENABLE_MPI==1
                _invmapper[gidMap].rank.push_back(m_rank);
#endif
            }
            indRef++;
            indMap++;
        }
        else if (octantsIRMapped->at(indMap).octant.getLevel() > octantsIRReference->at(indRef).octant.getLevel()){
            m_mapper[idRef].type = bitpit::mapping::Type::TYPE_COARSENING;

            uint64_t mortonlastdesc = meshReference->getTree().getLastDescMorton(&octantsIRReference->at(indRef).octant);
            uint64_t mortonMap = meshMapped->getTree().getMorton(&octantsIRMapped->at(indMap).octant);
            uint64_t mortonMapLastDesc = meshMapped->getTree().getLastDescMorton(&octantsIRMapped->at(indMap).octant);
            while(mortonMap <= mortonlastdesc){
                m_mapper[idRef].mapped.push_back(idMap);
#if BITPIT_ENABLE_MPI==1
                m_mapper[idRef].rank.push_back(rank);
#endif
                if (fillInv){
                    _invmapper[gidMap].type = bitpit::mapping::Type::TYPE_REFINEMENT;
                    _invmapper[gidMap].mapped.push_back(idRef);
#if BITPIT_ENABLE_MPI==1
                    _invmapper[gidMap].rank.push_back(m_rank);
#endif
                }
                indMap++;
                if (indMap == nMap) break;
                idMap = octantsIRMapped->at(indMap).ID;
                gidMap = octantsIRMapped->at(indMap).globalID;
#if BITPIT_ENABLE_MPI==1
                rank = octantsIRMapped->at(indMap).rank;
#endif
                mortonMap = meshMapped->getTree().getMorton(&octantsIRMapped->at(indMap).octant);
                mortonMapLastDesc = meshMapped->getTree().getLastDescMorton(&octantsIRMapped->at(indMap).octant);
            }
            if (mortonMapLastDesc >= mortonlastdesc) indRef++;
        }
        else if (octantsIRMapped->at(indMap).octant.getLevel() < octantsIRReference->at(indRef).octant.getLevel()){

            uint64_t morton = meshReference->getTree().getMorton(&octantsIRReference->at(indRef).octant);
            uint64_t mortonlastdescmesh = meshMapped->getTree().getLastDescMorton(&octantsIRMapped->at(indMap).octant);

            if (fillInv){
                _invmapper[gidMap].type = bitpit::mapping::Type::TYPE_COARSENING;
            }

            while (morton <= mortonlastdescmesh){
                m_mapper[idRef].type = bitpit::mapping::Type::TYPE_REFINEMENT;
                m_mapper[idRef].mapped.push_back(idMap);
#if BITPIT_ENABLE_MPI==1
                m_mapper[idRef].rank.push_back(rank);
#endif
                if (fillInv){
                    _invmapper[gidMap].mapped.push_back(idRef);
#if BITPIT_ENABLE_MPI==1
                    _invmapper[gidMap].rank.push_back(m_rank);
#endif
                }
                indRef++;
                if (indRef == nRef) break;
                morton = meshReference->getTree().getMorton(&octantsIRReference->at(indRef).octant);
                idRef = octantsIRReference->at(indRef).ID;
            }
            indMap++;
        }
    }

#if BITPIT_ENABLE_MPI==1
    if (meshMapped->isPartitioned() && !localMapped && !checkPartition())
        _communicateInverseMapper(_invmapper, octantsIRMapped);
    //After communication _invmapper is changed and it has info for each local ID as in serial
#endif

    for (pair<const long, mapping::MappingInfo> & val : _invmapper){
        m_invmapper[val.first] = val.second;
    }

}

#if BITPIT_ENABLE_MPI==1
/**
 * Initializes the MPI communicator to be used for parallel communications.
 *
 * \param communicator is the communicator.
 */
void MapperVolOctree::initializeCommunicator(MPI_Comm communicator)
{
    // Communication can be set just once
    if (isCommunicatorSet())
        throw std::runtime_error ("MapperVolOctree communicator can be set just once");

    // The communicator has to be valid
    if (communicator == MPI_COMM_NULL)
        throw std::runtime_error ("MapperVolOctree communicator is not valid");

    // Create a copy of the user-specified communicator
    //
    // No library routine should use MPI_COMM_WORLD as the communicator;
    // instead, a duplicate of a user-specified communicator should always
    // be used.
    MPI_Comm_dup(communicator, &m_communicator);
}

/**
 * Returns the MPI communicator stored within LevelSetKernel.
 * @return MPI communicator.
 */
MPI_Comm MapperVolOctree::getCommunicator() const
{
    return m_communicator;
}

/**
 * Checks if the communicator to be used for parallel communications has
 * already been set.
 *
 * \result Returns true if the communicator has been set, false otherwise.
 */
bool MapperVolOctree::isCommunicatorSet() const
{

    return (getCommunicator() != MPI_COMM_NULL);
}

/**
 * Frees the MPI communicator associated to the patch
 */
void MapperVolOctree::freeCommunicator()
{
    if (!isCommunicatorSet())
        return;

    int finalizedCalled;
    MPI_Finalized(&finalizedCalled);
    if (finalizedCalled)
        return;

    MPI_Comm_free(&m_communicator);
}

/**
 * Check if the reference and the mapped meshes have the same geometric partitioning
 *
 * \return True if the meshes have identical geometric partitioning, false otherwise
 *
 */
bool MapperVolOctree::checkPartition()
{
    bool check;

    std::vector<uint64_t> partitionReference = m_referenceMesh->getTree().getPartitionLastDesc();
    std::vector<uint64_t> partitionMapped = m_mappedMesh->getTree().getPartitionLastDesc();

    check = true;
    for (int irank = 0; irank < m_nProcs; irank++){
        check = check && (partitionReference[irank] == partitionMapped[irank]);
    }

    return check;
}

/**
 * Mapping computing between two partitioned meshes (stored in class members).
 *
 * \param[in] fillInv True if the inverse mapper has to be computed
 */
void MapperVolOctree::_mapMeshPartitioned(bool fillInv)
{
    _recoverPartition();
    long indRef = 0;
    //Fill IR with reference mesh (TODO make a method to do that)
    std::vector<OctantIR> octantsIRReference;
    long n = static_cast<VolOctree*>(m_referenceMesh)->getInternalCount();
    const std::vector<Octant>* octants = static_cast<VolOctree*>(m_referenceMesh)->getTree().getInternalOctants();
    for (long i = 0; i < n; i++){
        VolOctree::OctantInfo octinfoRef(i, true);
        long idRef = static_cast<VolOctree*>(m_referenceMesh)->getOctantId(octinfoRef);
        OctantIR val(octants->at(i), idRef, m_rank);
        octantsIRReference.push_back(val);
    }
    _mapMeshesSamePartition(&octantsIRReference, &m_partitionIR.list_rec_octantIR_before, fillInv, indRef);
    _mapMeshesSamePartition(&octantsIRReference, nullptr, fillInv, indRef);
    _mapMeshesSamePartition(&octantsIRReference, &m_partitionIR.list_rec_octantIR_after, fillInv, indRef);

}

/**
 * Recover the overlapped partitions between the mapped meshes.
 *
 * \return False if the partition structures are different from the ones stored in class (e.g. after an adaptation or after a load balancing).
 *
 */
bool MapperVolOctree::_recoverPartition()
{

    //Now the mapped mesh is repartitioned (is copied...) over the processes to match
    // the partitioning between the two meshes

    //find owner of reference partition over mapped mesh partition
    VolOctree* meshReference = static_cast<VolOctree*>(m_referenceMesh);
    VolOctree* meshMapped = static_cast<VolOctree*>(m_mappedMesh);
    std::vector<uint64_t> partitionLDReference = meshReference->getTree().getPartitionLastDesc();
    std::vector<uint64_t> partitionLDMapped = meshMapped->getTree().getPartitionLastDesc();
    std::vector<uint64_t> partitionFDReference = meshReference->getTree().getPartitionFirstDesc();
    std::vector<uint64_t> partitionFDMapped = meshMapped->getTree().getPartitionFirstDesc();

    if (m_partitionIR.partitionFDReference.size() != 0){
        if (m_partitionIR.partitionFDMapped != partitionFDMapped ||
                m_partitionIR.partitionLDMapped != partitionLDMapped ||
                m_partitionIR.partitionFDReference != partitionFDReference ||
                m_partitionIR.partitionLDReference != partitionLDReference){
            m_partitionIR.partitionFDReference.clear();
            m_partitionIR.partitionLDReference.clear();
            m_partitionIR.partitionFDMapped.clear();
            m_partitionIR.partitionLDMapped.clear();
            return false;
        }
    }

    std::map<int, std::vector<int>> frommapped_rank;
    std::map<int,std::vector<int>> toreference_rank;
    for (int ref_rank = 0; ref_rank < m_nProcs; ref_rank++){

        uint64_t local_first_morton = partitionFDReference[ref_rank];
        uint64_t local_last_morton = partitionLDReference[ref_rank];

        for (int irank = 0; irank < m_nProcs; irank++){
            if (irank != ref_rank){
                if (partitionFDMapped[irank] <= local_first_morton && partitionLDMapped[irank] >= local_first_morton){
                    frommapped_rank[ref_rank].push_back(irank);
                    toreference_rank[irank].push_back(ref_rank);
                }
                else if (partitionFDMapped[irank] <= local_last_morton && partitionLDMapped[irank] >= local_last_morton){
                    frommapped_rank[ref_rank].push_back(irank);
                    toreference_rank[irank].push_back(ref_rank);
                }
                else if (partitionFDMapped[irank] <= local_first_morton && partitionLDMapped[irank] >= local_last_morton){
                    frommapped_rank[ref_rank].push_back(irank);
                    toreference_rank[irank].push_back(ref_rank);
                }
                else if (partitionFDMapped[irank] >= local_first_morton && partitionLDMapped[irank] <= local_last_morton){
                    frommapped_rank[ref_rank].push_back(irank);
                    toreference_rank[irank].push_back(ref_rank);
                }
            }
        }
    }

    clearPartitionMappingLists();

    uint8_t dim = meshMapped->getDimension();

    //locally the mapped mesh build the lists of octants to send
    std::map<int, std::vector<Octant>> list_octant;
    std::map<int, std::vector<long>> list_ID;
    std::map<int, std::vector<long>> list_globalID;
    uint32_t idx = 0;
    uint64_t morton = meshMapped->getTree().getLastDescMorton(idx);
    for (int reference_rank : toreference_rank[m_rank]){
        uint64_t reference_first_morton = partitionFDReference[reference_rank];
        uint64_t reference_last_morton = partitionLDReference[reference_rank];
        while (morton < reference_first_morton){
            idx++;
            if (idx == meshMapped->getTree().getNumOctants()) break;
            morton = meshMapped->getTree().getLastDescMorton(idx);
        }
        while (morton < reference_last_morton){
            Octant oct = *meshMapped->getTree().getOctant(idx);
            list_octant[reference_rank].push_back(oct);
            VolOctree::OctantInfo octinfo(idx, true);
            long ID = meshMapped->getOctantId(octinfo);
            list_ID[reference_rank].push_back(ID);
            long globalID = meshMapped->getTree().getGlobalIdx(idx);
            list_globalID[reference_rank].push_back(globalID);
            OctantIR val(oct, ID, globalID, reference_rank);
            m_partitionIR.list_sent_octantIR.push_back(val);
            idx++;
            if (idx == meshMapped->getTree().getNumOctants()) break;
            morton = meshMapped->getTree().getMorton(idx);
        }
    }

    // communicate: send local mapped octants to reference rank and receive mapped octants
    //              from other processes

    //build send buffers
    DataCommunicator octCommunicator(m_communicator);
    MPI_Barrier(m_communicator);

    //set size
    //TODO make accessible global variables in ParaTree
    uint32_t x,y,z;
    long ID, globalID;
    uint8_t l;
    uint8_t octantBytes = uint8_t(sizeof(uint32_t)*3 + sizeof(uint8_t) + sizeof(long)*2);
    for (int reference_rank : toreference_rank[m_rank]){

        std::size_t buffSize = list_octant[reference_rank].size() * (std::size_t)octantBytes;
        octCommunicator.setSend(reference_rank,buffSize);

        //fill buffer with octants
        SendBuffer &sendBuffer = octCommunicator.getSendBuffer(reference_rank);
        int8_t m;
        bool info[17];
        uint32_t ii = 0;
        for(Octant & octant : list_octant[reference_rank]){
            x = octant.getX();
            y = octant.getY();
            z = octant.getZ();
            l = octant.getLevel();
            ID = list_ID[reference_rank][ii];
            globalID = list_globalID[reference_rank][ii];
            sendBuffer << x;
            sendBuffer << y;
            sendBuffer << z;
            sendBuffer << l;
            sendBuffer << ID;
            sendBuffer << globalID;
            ii++;
        }
    }

    octCommunicator.discoverRecvs();
    octCommunicator.startAllRecvs();
    octCommunicator.startAllSends();

    vector<int> recvRanks = octCommunicator.getRecvRanks();
    std::sort(recvRanks.begin(),recvRanks.end());

    std::vector<OctantIR> & list_rec_octantIR_before = m_partitionIR.list_rec_octantIR_before;
    std::vector<OctantIR> & list_rec_octantIR_after = m_partitionIR.list_rec_octantIR_after;
    std::vector<OctantIR>* _list_rec_octantIR;

    for(int rank : recvRanks){

        octCommunicator.waitRecv(rank);

        RecvBuffer & recvBuffer = octCommunicator.getRecvBuffer(rank);
        long bufferSize = recvBuffer.getSize();
        uint32_t nof = (uint32_t)(bufferSize / (uint32_t)octantBytes);
        if (rank < m_rank){
            _list_rec_octantIR = &list_rec_octantIR_before;
        }
        else if (rank > m_rank){
            _list_rec_octantIR = &list_rec_octantIR_after;
        }
        else
            break;

        for(int i = 0; i < nof; i++){
            recvBuffer >> x;
            recvBuffer >> y;
            recvBuffer >> z;
            recvBuffer >> l;
            recvBuffer >> ID;
            recvBuffer >> globalID;
            OctantIR val(Octant(dim,l,x,y,z), ID, globalID, rank);
            _list_rec_octantIR->push_back(val);
        }
        for (OctantIR & octir : *_list_rec_octantIR){
            ID = octir.ID;
            m_partitionIR.map_rank_rec_octantIR[rank][ID] = &octir;
        }
    }

    octCommunicator.waitAllSends();

    m_partitionIR.partitionLDReference = partitionLDReference;
    m_partitionIR.partitionLDMapped = partitionLDMapped;
    m_partitionIR.partitionFDReference = partitionFDReference;
    m_partitionIR.partitionFDMapped = partitionFDMapped;

    return true;
}

/**
 * Communicate inverse mapper info of overlapped partitions between processes.
 *
 * \param[in] _invmapper Temporary inverse mapper structure
 * \param[in] octantsIRMapped List of overlapped octants of mapped mesh
 *
 */
void MapperVolOctree::_communicateInverseMapper(std::unordered_map<long, mapping::MappingInfo> & _invmapper, const std::vector<OctantIR> * octantsIRMapped)
{

    //Recover mapper elements to send (to partitions of mapped mesh)
    std::set<int> toRanks;
    std::map<int, std::vector<long> > toRankgID;
    std::map<long, long > gIDtoID;
    for (const OctantIR & val : *octantsIRMapped)
    {
        toRankgID[val.rank].push_back(val.globalID);
        gIDtoID[val.globalID] = val.ID;
        toRanks.insert(val.rank);
    }

    //build send buffers
    DataCommunicator mapCommunicator(m_communicator);
    MPI_Barrier(m_communicator);

    //set size
    for (int rank : toRanks){
        std::size_t buffSize = 0;
        int mappedSize;
        for (long & gID : toRankgID[rank]){
            mapping::MappingInfo _info = _invmapper[gID];
            mappedSize = _info.mapped.size();
            std::size_t infoBytes = std::size_t(sizeof(long) + sizeof(int) + sizeof(int) + sizeof(int) + (sizeof(long))*mappedSize);
            buffSize += infoBytes;
        }
        buffSize += std::size_t(sizeof(int));
        mapCommunicator.setSend(rank,buffSize);

        //fill buffer with octants and local map for inverse mapper in PartitionMapper member
        SendBuffer &sendBuffer = mapCommunicator.getSendBuffer(rank);
        mapping::Type type;
        mapping::Entity entity;
        std::vector<long> mapped;
        sendBuffer << int(toRankgID[rank].size());
        for (long & gID : toRankgID[rank]){
            mapping::MappingInfo _info = _invmapper[gID];
            sendBuffer << gID;
            sendBuffer << int(_info.type);
            sendBuffer << int(_info.entity);
            int mappedSize = _info.mapped.size();
            sendBuffer << mappedSize;
            for (long refId : _info.mapped){
                sendBuffer << refId;
            }

            m_partitionIR.map_rank_invmapper[rank][gIDtoID[gID]] = _info;

        }
    }

    mapCommunicator.discoverRecvs();
    mapCommunicator.startAllRecvs();
    mapCommunicator.startAllSends();

    vector<int> recvRanks = mapCommunicator.getRecvRanks();
    std::sort(recvRanks.begin(),recvRanks.end());

    _invmapper.clear();

    for(int rank : recvRanks){

        mapCommunicator.waitRecv(rank);

        RecvBuffer & recvBuffer = mapCommunicator.getRecvBuffer(rank);
        long bufferSize = recvBuffer.getSize();

        int nof;
        recvBuffer >> nof;

        for(std::size_t i = 0; i < nof; i++){
            mapping::MappingInfo _info;
            long localID;
            long globalID;
            int type;
            int entity;
            recvBuffer >> globalID;
            uint32_t idx = static_cast<VolOctree*>(m_mappedMesh)->getTree().getLocalIdx(globalID);
            VolOctree::OctantInfo octinfo(idx, true);
            localID = static_cast<VolOctree*>(m_mappedMesh)->getOctantId(octinfo);

            recvBuffer >> type;
            recvBuffer >> entity;
            _info.type = mapping::Type(type);
            _info.entity = mapping::Entity(entity);
            int nmap;
            recvBuffer >> nmap;
            for (std::size_t j=0; j<nmap; j++){
                long idref;
                recvBuffer >> idref;
                _info.mapped.push_back(idref);
                _info.rank.push_back(rank);
            }
            _invmapper[localID] = _info;
        }
    }

    mapCommunicator.waitAllSends();

}

/**
 * Communicate inverse mapper info of overlapped partitions back to the processes of mapped partitions.
 *
 */
void MapperVolOctree::_communicateInverseMapperBack()
{

    //Recover mapper elements to send (to partitions of mapped mesh)
    std::set<int> toRanks;
    std::map<int, std::vector<long> > toRankID;
    for (const std::pair<int, std::unordered_map<long, mapping::MappingInfo> > & val : m_partitionIR.map_rank_invmapper)
    {
        toRanks.insert(val.first);

        std::unordered_map<long, mapping::MappingInfo> _invmapper = m_partitionIR.map_rank_invmapper[val.first];
        for (const std::pair<long, mapping::MappingInfo> & sub : _invmapper)
            toRankID[val.first].push_back(sub.first);
    }

    //build send buffers
    DataCommunicator mapCommunicator(m_communicator);
    MPI_Barrier(m_communicator);

    //set size
    for (int rank : toRanks){
        std::size_t buffSize = 0;
        int mappedSize;
        std::unordered_map<long, mapping::MappingInfo> _invmapper = m_partitionIR.map_rank_invmapper[rank];
        for (long & ID : toRankID[rank]){
            mapping::MappingInfo _info = _invmapper[ID];
            mappedSize = _info.mapped.size();
            std::size_t infoBytes = std::size_t(sizeof(long) + sizeof(int) + sizeof(int) + sizeof(int) + (sizeof(long))*mappedSize);
            buffSize += infoBytes;
        }
        buffSize += std::size_t(sizeof(int));
        mapCommunicator.setSend(rank,buffSize);

        //fill buffer with octants and local map for inverse mapper in PartitionMapper member
        SendBuffer &sendBuffer = mapCommunicator.getSendBuffer(rank);
        std::vector<long> mapped;
        sendBuffer << int(toRankID[rank].size());
        for (long & ID : toRankID[rank]){
            mapping::MappingInfo _info = _invmapper[ID];
            sendBuffer << ID;
            sendBuffer << int(_info.type);
            sendBuffer << int(_info.entity);
            int mappedSize = _info.mapped.size();
            sendBuffer << mappedSize;
            for (long refId : _info.mapped){
                sendBuffer << refId;
            }
        }
    }

    MPI_Barrier(m_communicator);

    mapCommunicator.discoverRecvs();
    mapCommunicator.startAllRecvs();
    mapCommunicator.startAllSends();

    vector<int> recvRanks = mapCommunicator.getRecvRanks();
    std::sort(recvRanks.begin(),recvRanks.end());

    for(int rank : recvRanks){

        mapCommunicator.waitRecv(rank);

        RecvBuffer & recvBuffer = mapCommunicator.getRecvBuffer(rank);
        int nof;
        recvBuffer >> nof;
        for(std::size_t i = 0; i < nof; i++){
            mapping::MappingInfo _info;
            long localID;
            int type;
            int entity;
            recvBuffer >> localID;
            recvBuffer >> type;
            recvBuffer >> entity;
            _info.type = mapping::Type(type);
            _info.entity = mapping::Entity(entity);
            int nmap;
            recvBuffer >> nmap;
            for (std::size_t j=0; j<nmap; j++){
                long idref;
                recvBuffer >> idref;
                _info.mapped.push_back(idref);
                _info.rank.push_back(rank);
            }
            m_invmapper[localID] = _info;
        }
    }

    mapCommunicator.waitAllSends();

}

/**
 * Communicate adaption info of overlapped partitions of the mapped mesh to the processes of reference partitions.
 *
 * \param[in] infoAdapt Adaption info of local mapped partitions
 * \param[out] MappingInfoAdapt Adaption info to local reference partitions
 *
 */
void MapperVolOctree::_communicateMappedAdaptionInfo(const std::vector<adaption::Info> & infoAdapt, std::vector<adaption::Info> & MappingInfoAdapt)
{

    //Recover mapper elements to send (to partitions of mapped mesh)
    std::set<int> toRanks;
    std::map<int, std::vector<long> > toRankID;
    std::map<int, std::set<long> > toRankInd;
    int ii = 0;
    for (auto & info : infoAdapt){
        for (const long & ID : info.previous){
            for (const int & rank : m_previousmapper[ID].rank){
                toRanks.insert(rank);
                toRankInd[rank].insert(ii);
                toRankID[rank].push_back(ID);
            }
        }
        ii++;
    }

    //build send buffers
    DataCommunicator mapCommunicator(m_communicator);
    MPI_Barrier(m_communicator);

    //set size previousmapper + infoadapt
    for (int rank : toRanks){
        std::size_t buffSize = 0;
        int mappedSize;
        for (long & ID : toRankID[rank]){
            mapping::MappingInfo _info = m_previousmapper[ID];
            mappedSize = _info.mapped.size();
            std::size_t infoBytes = std::size_t(sizeof(long) + 3*sizeof(int) + (sizeof(long))*mappedSize);
            buffSize += infoBytes;
        }
        buffSize += std::size_t(sizeof(int));

        int currentSize, previousSize;
        for (const long & ind : toRankInd[rank]){
            adaption::Info _info = infoAdapt[ind];
            currentSize = _info.current.size();
            previousSize = _info.previous.size();
            std::size_t infoBytes = std::size_t(5*sizeof(int) + (sizeof(long))*(currentSize+previousSize));
            buffSize += infoBytes;
        }
        buffSize += std::size_t(sizeof(int));

        mapCommunicator.setSend(rank,buffSize);

        //fill buffer with octants and local map for inverse mapper in PartitionMapper member
        SendBuffer &sendBuffer = mapCommunicator.getSendBuffer(rank);
        std::vector<long> mapped;
        sendBuffer << int(toRankID[rank].size());
        for (long & ID : toRankID[rank]){
            mapping::MappingInfo _info = m_previousmapper[ID];
            sendBuffer << ID;
            sendBuffer << int(_info.type);
            sendBuffer << int(_info.entity);
            int mappedSize = _info.mapped.size();
            sendBuffer << mappedSize;
            for (long refId : _info.mapped){
                sendBuffer << refId;
            }
        }

        std::vector<long> current;
        std::vector<long> previous;
        sendBuffer << int(toRankInd[rank].size());
        for (const long & ind : toRankInd[rank]){
            adaption::Info _info = infoAdapt[ind];
            sendBuffer << int(_info.type);
            sendBuffer << int(_info.entity);
            int currentSize = _info.current.size();
            int previousSize = _info.previous.size();
            sendBuffer << currentSize;
            for (long Id : _info.current){
                sendBuffer << Id;
            }
            sendBuffer << previousSize;
            for (long Id : _info.previous){
                sendBuffer << Id;
            }
            //rank of adaption is used to identify the rank where the adaption occurs
            sendBuffer << m_rank;
        }

    }

    mapCommunicator.discoverRecvs();
    mapCommunicator.startAllRecvs();
    mapCommunicator.startAllSends();

    vector<int> recvRanks = mapCommunicator.getRecvRanks();
    std::sort(recvRanks.begin(),recvRanks.end());

    for(int rank : recvRanks){

        mapCommunicator.waitRecv(rank);

        RecvBuffer & recvBuffer = mapCommunicator.getRecvBuffer(rank);
        int nofMapper;
        recvBuffer >> nofMapper;
        for(std::size_t i = 0; i < nofMapper; i++){
            mapping::MappingInfo _info;
            long localID;
            int mtype;
            int mentity;
            recvBuffer >> localID;
            recvBuffer >> mtype;
            recvBuffer >> mentity;
            _info.type = mapping::Type(mtype);
            _info.entity = mapping::Entity(mentity);
            int nmap;
            recvBuffer >> nmap;
            for (std::size_t j=0; j<nmap; j++){
                long idref;
                recvBuffer >> idref;
                _info.mapped.push_back(idref);
                _info.rank.push_back(m_rank);
            }
            m_partitionIR.map_rank_previousmapper[rank][localID] = _info;
            m_partitionIR.map_rank_invmapper[rank].erase(localID);
        }

        int nofInfo;
        recvBuffer >> nofInfo;
        for(std::size_t i = 0; i < nofInfo; i++){
            adaption::Info _info;
            int type;
            int entity;
            recvBuffer >> type;
            recvBuffer >> entity;
            _info.type = adaption::Type(type);
            _info.entity = adaption::Entity(entity);
            int ncurrent;
            recvBuffer >> ncurrent;
            for (std::size_t j=0; j<ncurrent; j++){
                long id;
                recvBuffer >> id;
                _info.current.push_back(id);
            }
            int nprevious;
            recvBuffer >> nprevious;
            for (std::size_t j=0; j<nprevious; j++){
                long id;
                recvBuffer >> id;
                _info.previous.push_back(id);
            }
            int arank;
            recvBuffer >> arank;
            _info.rank = arank;
            MappingInfoAdapt.push_back(_info);
        }

    }

    mapCommunicator.waitAllSends();

}

#endif

}
