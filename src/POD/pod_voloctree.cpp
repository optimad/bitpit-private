/*---------------------------------------------------------------------------*\
 *
 *  bitpit
 *
 *  Copyright (C) 2015-2019 OPTIMAD engineering Srl
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

#include "pod_voloctree.hpp"

namespace bitpit {

/**
 * \class PODVolOctree
 * \ingroup POD
 *
 * \brief The PODVolOctree is the specialized class of PODKernel for VolOctree meshes.
 *
 * The PODVolOctree is the specialized class of PODKernel for VolOctree meshes.
 */

/**
 * Creates a new PODVolOctree object.
 */
# if BITPIT_ENABLE_MPI
/**
 * \param[in] comm The MPI communicator used by the pod object. MPI_COMM_WORLD is the default value.
 */
PODVolOctree::PODVolOctree(MPI_Comm comm):PODKernel(comm)
# else
PODVolOctree::PODVolOctree():PODKernel()
# endif
{
}

/**
 * Destructor of PODVolOctree
 */
PODVolOctree::~PODVolOctree()
{
}

/**
 * Create the mesh for POD.
 *
 * \return Pointer to instantiated mesh.
 */
VolumeKernel* PODVolOctree::createMesh()
{
    VolumeKernel *mesh = new VolOctree();

#if BITPIT_ENABLE_MPI
    mesh->setCommunicator(getCommunicator());
#endif

    return mesh;
}

/**
 * Compute the mapping of an input mesh on the pod mesh.
 * \param[in] mesh Pointer to input mesh.
 * \param[in] fillInv If true even the inverse mapping is computed.
 */
void PODVolOctree::_computeMapper(VolumeKernel* mesh, MapperVolOctree* & mapper, bool fillInv)
{
    VolOctree* meshPOD = static_cast<VolOctree*>(getMesh());
    VolOctree* _mesh = static_cast<VolOctree*>(mesh);

    mapper = new MapperVolOctree(meshPOD, _mesh, getCommunicator());
    mapper->initialize(fillInv);
}

/**
 * Map a POD Field TO the POD mesh.
 *
 * \param[in] field Input POD Field
 * \param[in] targetCells Pointer to set of target cells of POD mesh (cells with mapped values as result)
 * \return Mapped POD field on target cells of POD mesh.
 */
pod::PODField PODVolOctree::mapPODFieldToPOD(const pod::PODField & field, const std::unordered_set<long> * targetCells)
{

    //Check target cells
    std::unordered_set<long> targetCellsStorage;
    if (!targetCells) {
        for (const Cell &cell : getMesh()->getCells())
            targetCellsStorage.insert(cell.getId());

        targetCells = &targetCellsStorage;
    }

    // Map data fields on pod mesh
    std::size_t nsf = field.scalar->getFieldCount();
    std::size_t nvf = field.vector->getFieldCount();

    const PiercedStorage<mapping::Info> & m_mapper = getMapper()->getMapping();

    pod::PODField mappedField(nsf, nvf, getMesh(), &getMesh()->getCells());


#if BITPIT_ENABLE_MPI
    if ( getMapper()->checkPartition()){
#endif

        for (long id : *targetCells){
            std::size_t rawIndex = m_meshPOD->getCells().getRawIndex(id);

            double *datamappedS = mappedField.scalar->data(id);
            std::array<double,3> *datamappedV = mappedField.vector->data(id);

            if (m_mapper[id].type == mapping::Type::TYPE_RENUMBERING){
                bool dataB = field.mask->at(m_mapper[id].mapped[0]);
                mappedField.mask->set(id, dataB);

                double *dataS = field.scalar->data(m_mapper[id].mapped[0]);
                for (std::size_t i = 0; i < nsf; i++){
                    double *dataSi = dataS + i;
                    double *datamappedSi = datamappedS + i;
                    (*datamappedSi) = (*dataSi);
                }

                std::array<double,3> *dataV = field.vector->data(m_mapper[id].mapped[0]);
                for (std::size_t i = 0; i < nvf; i++) {
                    std::array<double,3> *dataVi = dataV + i;
                    std::array<double,3> *datamappedVi = datamappedV + i;
                    (*datamappedVi) = (*dataVi);
                }
            }
            else if (m_mapper[id].type == mapping::Type::TYPE_COARSENING){
                mappedField.mask->set(id, true);

                for (std::size_t i = 0; i < nsf; i++){
                    double *datamappedSi = datamappedS + i;
                    (*datamappedSi) = 0.0;
                }
                for (std::size_t i = 0; i < nvf; i++){
                    std::array<double,3> *datamappedVi = datamappedV + i;
                    (*datamappedVi) = std::array<double,3>{0.0, 0.0, 0.0};
                }

                bool dataB, dataMappedB = true;
                double *dataS;
                std::array<double,3> *dataV;
                double volmapped = getRawCellVolume(rawIndex);
                for (long idd : m_mapper[id].mapped){
                    dataB = field.mask->at(idd);
                    dataMappedB &= dataB;

                    dataS = field.scalar->data(idd);
                    double vol = field.mesh->evalCellVolume(idd);
                    for (std::size_t i = 0; i < nsf; i++){
                        double *dataSi = dataS + i;
                        double *datamappedSi = datamappedS + i;
                        (*datamappedSi) += (*dataSi) * vol / volmapped;
                    }
                    dataV = field.vector->data(idd);
                    for (std::size_t i = 0; i < nvf; i++) {
                        std::array<double,3> *dataVi = dataV + i;
                        std::array<double,3> *datamappedVi = datamappedV + i;
                        (*datamappedVi) += (*dataVi) * vol / volmapped;
                    }
                }
                mappedField.mask->set(id, dataMappedB);

            }
            else if (m_mapper[id].type == mapping::Type::TYPE_REFINEMENT){
                bool dataB = field.mask->at(m_mapper[id].mapped[0]);
                double *dataS = field.scalar->data(m_mapper[id].mapped[0]);
                std::array<double,3> *dataV = field.vector->data(m_mapper[id].mapped[0]);

                mappedField.mask->set(id, dataB);

                for (std::size_t i = 0; i < nsf; i++){
                    double *dataSi = dataS + i;
                    double *datamappedSi = datamappedS + i;
                    (*datamappedSi) = (*dataSi);
                }
                for (std::size_t i = 0; i < nvf; i++) {
                    std::array<double,3> *dataVi = dataV + i;
                    std::array<double,3> *datamappedVi = datamappedV + i;
                    (*datamappedVi) = (*dataVi);
                }
            }
        }

#if BITPIT_ENABLE_MPI
    }
    else{

        //Mapping between partitioned meshes

        //Communicate data and volumes
        std::map<int, std::map<long, bool> > dataBrec;
        std::map<int, std::map<long, std::vector<double> > > dataSrec;
        std::map<int, std::map<long, std::vector<std::array<double, 3> > > > dataVrec;
        std::map<int, std::map<long, double> > volrec;

        communicatePODField(field, dataBrec, dataSrec, dataVrec, volrec);

        for (long id : *targetCells){

            double *datamappedS = mappedField.scalar->data(id);
            std::array<double,3> *datamappedV = mappedField.vector->data(id);

            if (m_mapper[id].type == mapping::Type::TYPE_RENUMBERING){

                int rank = m_mapper[id].rank[0];
                if (rank != m_rank){
                    bool dataB = dataBrec[rank][m_mapper[id].mapped[0]];
                    mappedField.mask->set(id, dataB);

                    double *dataS = dataSrec[rank][m_mapper[id].mapped[0]].data();
                    for (std::size_t i = 0; i < nsf; i++){
                        double *dataSi = dataS + i;
                        double *datamappedSi = datamappedS + i;
                        (*datamappedSi) = (*dataSi);
                    }

                    std::array<double,3> *dataV = dataVrec[rank][m_mapper[id].mapped[0]].data();
                    for (std::size_t i = 0; i < nvf; i++) {
                        std::array<double,3> *dataVi = dataV + i;
                        std::array<double,3> *datamappedVi = datamappedV + i;
                        (*datamappedVi) = (*dataVi);
                    }

                }
                else{
                    bool dataB = field.mask->at(m_mapper[id].mapped[0]);
                    mappedField.mask->set(id, dataB);

                    double *dataS = field.scalar->data(m_mapper[id].mapped[0]);
                    for (std::size_t i = 0; i < nsf; i++){
                        double *dataSi = dataS + i;
                        double *datamappedSi = datamappedS + i;
                        (*datamappedSi) = (*dataSi);
                    }

                    std::array<double,3> *dataV = field.vector->data(m_mapper[id].mapped[0]);
                    for (std::size_t i = 0; i < nvf; i++) {
                        std::array<double,3> *dataVi = dataV + i;
                        std::array<double,3> *datamappedVi = datamappedV + i;
                        (*datamappedVi) = (*dataVi);
                    }
                }//end if other rank

            }
            else if (m_mapper[id].type == mapping::Type::TYPE_COARSENING){
                mappedField.mask->set(id, true);

                for (std::size_t i = 0; i < nsf; i++){
                    double *datamappedSi = datamappedS + i;
                    (*datamappedSi) = 0.0;
                }
                for (std::size_t i = 0; i < nvf; i++){
                    std::array<double,3> *datamappedVi = datamappedV + i;
                    (*datamappedVi) = std::array<double,3>{0.0, 0.0, 0.0};
                }

                bool dataB, dataMappedB = false;
                double *dataS;
                std::array<double,3> *dataV;
                double vol;
                double volmapped = getCellVolume(id);
                int ii = 0;
                for (long idd : m_mapper[id].mapped){
                    int rank = m_mapper[id].rank[ii];
                    if (rank != m_rank){
                        dataB =  dataBrec[rank][idd];
                        dataS = dataSrec[rank][idd].data();
                        dataV = dataVrec[rank][idd].data();
                        vol = volrec[rank][idd];
                    }
                    else{
                        dataB =  field.mask->at(idd);
                        dataS = field.scalar->data(idd);
                        dataV = field.vector->data(idd);
                        vol = field.mesh->evalCellVolume(idd);
                    }
                    dataMappedB &= dataB;

                    for (std::size_t i = 0; i < nsf; i++){
                        double *dataSi = dataS + i;
                        double *datamappedSi = datamappedS + i;
                        (*datamappedSi) += (*dataSi) * vol / volmapped;
                    }
                    for (std::size_t i = 0; i < nvf; i++) {
                        std::array<double,3> *dataVi = dataV + i;
                        std::array<double,3> *datamappedVi = datamappedV + i;
                        (*datamappedVi) += (*dataVi) * vol / volmapped;
                    }
                    ii++;
                }
                mappedField.mask->set(id, dataMappedB);

            }
            else if (m_mapper[id].type == mapping::Type::TYPE_REFINEMENT){
                bool dataB;
                double *dataS;
                std::array<double,3> *dataV;
                int rank = m_mapper[id].rank[0];
                if (rank != m_rank){
                    dataB =  dataBrec[rank][m_mapper[id].mapped[0]];
                    dataS = dataSrec[rank][m_mapper[id].mapped[0]].data();
                    dataV = dataVrec[rank][m_mapper[id].mapped[0]].data();
                }
                else{
                    dataB =  field.mask->at(m_mapper[id].mapped[0]);
                    dataS = field.scalar->data(m_mapper[id].mapped[0]);
                    dataV = field.vector->data(m_mapper[id].mapped[0]);
                }

                mappedField.mask->set(id, dataB);

                for (std::size_t i = 0; i < nsf; i++){
                    double *dataSi = dataS + i;
                    double *datamappedSi = datamappedS + i;
                    (*datamappedSi) = (*dataSi);
                }
                for (std::size_t i = 0; i < nvf; i++) {
                    std::array<double,3> *dataVi = dataV + i;
                    std::array<double,3> *datamappedVi = datamappedV + i;
                    (*datamappedVi) = (*dataVi);
                }
            }
        } // end target cells

    }//end partitioned meshes
#endif
    return mappedField;

}

/**
 * Map a POD Field FROM the POD mesh.
 *
 * \param[out] field Output POD Field
 * \param[in] targetCells Pointer to set of target cells of output field mesh (cells with mapped values as result)
 * \param[in] mappedField POD field on target cells of output field mesh.
 */
void PODVolOctree::mapPODFieldFromPOD(pod::PODField & field, const std::unordered_set<long> * targetCells,
        const pod::PODField & mappedField)
{

    //Check target cells
    std::unordered_set<long> targetCellsStorage;
    if (!targetCells) {
        for (const Cell &cell : field.mesh->getCells())
            targetCellsStorage.insert(cell.getId());

        targetCells = &targetCellsStorage;
    }

    // Map data fields from pod mesh
    std::size_t nsf = field.scalar->getFieldCount();
    std::size_t nvf = field.vector->getFieldCount();

    const PiercedStorage<mapping::Info> & m_invmapper = getMapper()->getInverseMapping();

#if BITPIT_ENABLE_MPI
    if ( getMapper()->checkPartition()){
#endif


        for (long id : *targetCells){
            double *dataS = field.scalar->data(id);
            std::array<double,3> *dataV = field.vector->data(id);

            if (m_invmapper[id].type == mapping::Type::TYPE_RENUMBERING){
                bool datamappedB = mappedField.mask->at(m_invmapper[id].mapped[0]);
                field.mask->set(id, datamappedB);
                double *datamappedS = mappedField.scalar->data(m_invmapper[id].mapped[0]);
                for (std::size_t i = 0; i < nsf; i++){
                    double *dataSi = dataS + i;
                    double *datamappedSi = datamappedS + i;
                    (*dataSi) = (*datamappedSi);
                }
                std::array<double,3>  *datamappedV = mappedField.vector->data(m_invmapper[id].mapped[0]);
                for (std::size_t i = 0; i < nvf; i++) {
                    std::array<double,3> *dataVi = dataV + i;
                    std::array<double,3> *datamappedVi = datamappedV + i;
                    (*dataVi) = (*datamappedVi);
                }
            }
            else if (m_invmapper[id].type == mapping::Type::TYPE_COARSENING){
                bool dataB = true;
                for (std::size_t i = 0; i < nsf; i++){
                    double *dataSi = dataS + i;
                    (*dataSi) = 0.0;
                }
                for (std::size_t i = 0; i < nvf; i++){
                    std::array<double,3> *dataVi = dataV + i;
                    (*dataVi) = std::array<double,3>{0.0, 0.0, 0.0};
                }
                bool datamappedB;
                double *datamappedS;
                std::array<double,3> *datamappedV;
                double vol = field.mesh->evalCellVolume(id);
                for (long idd : m_invmapper[id].mapped){
                    std::size_t rawIndexIdd = m_meshPOD->getCells().getRawIndex(idd);
                    datamappedB = mappedField.mask->at(idd);
                    dataB &= datamappedB;
                    datamappedS = mappedField.scalar->data(idd);
                    double volmapped = getRawCellVolume(rawIndexIdd);
                    for (std::size_t i = 0; i < nsf; i++){
                        double *dataSi = dataS + i;
                        double *datamappedSi = datamappedS + i;
                        (*dataSi) += (*datamappedSi) * volmapped / vol;
                    }
                    datamappedV = mappedField.vector->data(idd);
                    for (std::size_t i = 0; i < nvf; i++) {
                        std::array<double,3> *dataVi = dataV + i;
                        std::array<double,3> *datamappedVi = datamappedV + i;
                        (*dataVi) += (*datamappedVi) * volmapped / vol;
                    }
                }
                field.mask->set(id, dataB);
            }
            else if (m_invmapper[id].type == mapping::Type::TYPE_REFINEMENT){
                bool datamappedB = mappedField.mask->at(m_invmapper[id].mapped[0]);
                field.mask->set(id, datamappedB);
                double *datamappedS = mappedField.scalar->data(m_invmapper[id].mapped[0]);
                for (std::size_t i = 0; i < nsf; i++){
                    double *dataSi = dataS + i;
                    double *datamappedSi = datamappedS + i;
                    (*dataSi) = (*datamappedSi);
                }
                std::array<double,3> *datamappedV = mappedField.vector->data(m_invmapper[id].mapped[0]);
                for (std::size_t i = 0; i < nvf; i++) {
                    std::array<double,3> *dataVi = dataV + i;
                    std::array<double,3> *datamappedVi = datamappedV + i;
                    (*dataVi) = (*datamappedVi);
                }
            }
        }//end target cells

#if BITPIT_ENABLE_MPI
    }
    else{

        //Mapping between partitioned meshes

        //Communicate data and volumes
        std::map<int, std::map<long, bool> > dataBrec;
        std::map<int, std::map<long, std::vector<double> > > dataSrec;
        std::map<int, std::map<long, std::vector<std::array<double, 3> > > > dataVrec;
        std::map<int, std::map<long, double> > volrec;

        communicatePODFieldFromPOD(mappedField, dataBrec, dataSrec, dataVrec, volrec);

        for (long id : *targetCells){
            double *dataS = field.scalar->data(id);
            std::array<double,3> *dataV = field.vector->data(id);

            if (m_invmapper[id].type == mapping::Type::TYPE_RENUMBERING){
                bool datamappedB;
                double *datamappedS;
                std::array<double,3>  *datamappedV;
                int rank = m_invmapper[id].rank[0];
                if (rank != m_rank){
                    datamappedB = dataBrec[rank][m_invmapper[id].mapped[0]];
                    datamappedS = dataSrec[rank][m_invmapper[id].mapped[0]].data();
                    datamappedV = dataVrec[rank][m_invmapper[id].mapped[0]].data();
                }
                else{
                    datamappedB = mappedField.mask->at(m_invmapper[id].mapped[0]);
                    datamappedS = mappedField.scalar->data(m_invmapper[id].mapped[0]);
                    datamappedV = mappedField.vector->data(m_invmapper[id].mapped[0]);
                }//endif different rank

                field.mask->set(id, datamappedB);
                for (std::size_t i = 0; i < nsf; i++){
                    double *dataSi = dataS + i;
                    double *datamappedSi = datamappedS + i;
                    (*dataSi) = (*datamappedSi);
                }
                for (std::size_t i = 0; i < nvf; i++) {
                    std::array<double,3> *dataVi = dataV + i;
                    std::array<double,3> *datamappedVi = datamappedV + i;
                    (*dataVi) = (*datamappedVi);
                }
            }
            else if (m_invmapper[id].type == mapping::Type::TYPE_COARSENING){
                bool dataB = true;
                for (std::size_t i = 0; i < nsf; i++){
                    double *dataSi = dataS + i;
                    (*dataSi) = 0.0;
                }
                for (std::size_t i = 0; i < nvf; i++){
                    std::array<double,3> *dataVi = dataV + i;
                    (*dataVi) = std::array<double,3>{0.0, 0.0, 0.0};
                }
                bool datamappedB;
                double *datamappedS;
                std::array<double,3> *datamappedV;
                double volmapped;
                double vol = field.mesh->evalCellVolume(id);
                int ii=0;
                for (long idd : m_invmapper[id].mapped){
                    int rank = m_invmapper[id].rank[ii];
                    if (rank != m_rank){
                        datamappedB = dataBrec[rank][idd];
                        datamappedS = dataSrec[rank][idd].data();
                        datamappedV = dataVrec[rank][idd].data();
                        volmapped = volrec[rank][idd];
                    }
                    else{
                        datamappedB = mappedField.mask->at(idd);
                        datamappedS = mappedField.scalar->data(idd);
                        datamappedV = mappedField.vector->data(idd);
                        volmapped = getCellVolume(idd);
                    }//endif different rank

                    dataB &= datamappedB;
                    for (std::size_t i = 0; i < nsf; i++){
                        double *dataSi = dataS + i;
                        double *datamappedSi = datamappedS + i;
                        (*dataSi) += (*datamappedSi) * volmapped / vol;
                    }
                    for (std::size_t i = 0; i < nvf; i++) {
                        std::array<double,3> *dataVi = dataV + i;
                        std::array<double,3> *datamappedVi = datamappedV + i;
                        (*dataVi) += (*datamappedVi) * volmapped / vol;
                    }
                    ii++;
                }
                field.mask->set(id, dataB);
            }
            else if (m_invmapper[id].type == mapping::Type::TYPE_REFINEMENT){
                bool datamappedB;
                double *datamappedS;
                std::array<double,3> *datamappedV;
                int rank = m_invmapper[id].rank[0];
                if (rank != m_rank){
                    datamappedB = dataBrec[rank][m_invmapper[id].mapped[0]];
                    datamappedS = dataSrec[rank][m_invmapper[id].mapped[0]].data();
                    datamappedV = dataVrec[rank][m_invmapper[id].mapped[0]].data();
                }
                else{
                    datamappedB = mappedField.mask->at(m_invmapper[id].mapped[0]);
                    datamappedS = mappedField.scalar->data(m_invmapper[id].mapped[0]);
                    datamappedV = mappedField.vector->data(m_invmapper[id].mapped[0]);
                }//endif different rank

                field.mask->set(id, datamappedB);
                for (std::size_t i = 0; i < nsf; i++){
                    double *dataSi = dataS + i;
                    double *datamappedSi = datamappedS + i;
                    (*dataSi) = (*datamappedSi);
                }
                for (std::size_t i = 0; i < nvf; i++) {
                    std::array<double,3> *dataVi = dataV + i;
                    std::array<double,3> *datamappedVi = datamappedV + i;
                    (*dataVi) = (*datamappedVi);
                }
            }
        }//end target cells

    }//endif partitioned
#endif

}

/**
 * Map a storage of fields TO the POD mesh.
 *
 * \param[in] field Input Fields storage
 * \param[in] mesh Pointer to input field mesh
 * \param[in] targetCells Pointer to set of target cells of POD mesh (cells with mapped values as result)
 * \param[in] scalarIds Ids of scalar fields in PiercedStorage.
 * \param[in] vectorIds Ids of vector fields in PiercedStorage.
 * \return Fields storage on target cells of POD mesh.
 */
PiercedStorage<double> PODVolOctree::mapFieldsToPOD(const PiercedStorage<double> & fields, const VolumeKernel * mesh,
        const std::unordered_set<long> * targetCells,
        const std::vector<std::size_t> &scalarIds,
        const std::vector<std::array<std::size_t, 3>> &vectorIds)
{

    //Check target cells
    std::unordered_set<long> targetCellsStorage;
    if (!targetCells) {
        for (const Cell &cell : getMesh()->getCells())
            targetCellsStorage.insert(cell.getId());

        targetCells = &targetCellsStorage;
    }

    // Map data fields on pod mesh
    std::size_t nsf = scalarIds.size();
    std::size_t nvf = vectorIds.size();

    const PiercedStorage<mapping::Info> & m_mapper = getMapper()->getMapping();

    PiercedStorage<double> mappedFields(fields.getFieldCount(), &getMesh()->getCells());

#if BITPIT_ENABLE_MPI
    if ( getMapper()->checkPartition()){
#endif

        for (long id : *targetCells){
            std::size_t rawIndex = m_meshPOD->getCells().getRawIndex(id);
            double *datamapped = mappedFields.data(id);
            if (m_mapper[id].type == mapping::Type::TYPE_RENUMBERING){
                const double *data = fields.data(m_mapper[id].mapped[0]);
                for (std::size_t i = 0; i < nsf; i++){
                    const double *datai = data + scalarIds[i];
                    double *datamappedi = datamapped + scalarIds[i];
                    (*datamappedi) = (*datai);
                }
                for (std::size_t i = 0; i < nvf; i++) {
                    for (std::size_t j = 0; j < 3; j++) {
                        const double *datai = data + vectorIds[i][j];
                        double *datamappedi = datamapped + vectorIds[i][j];
                        (*datamappedi) = (*datai);
                    }
                }
            }
            else if (m_mapper[id].type == mapping::Type::TYPE_COARSENING){
                for (std::size_t i = 0; i < nsf; i++){
                    double *datamappedi = datamapped + scalarIds[i];
                    (*datamappedi) = 0.0;
                }
                for (std::size_t i = 0; i < nvf; i++){
                    for (std::size_t j = 0; j < 3; j++) {
                        double *datamappedi = datamapped + vectorIds[i][j];
                        (*datamappedi) = 0.0;
                    }
                }
                double volmapped = getRawCellVolume(rawIndex);
                for (long idd : m_mapper[id].mapped){
                    const double *data = fields.data(idd);
                    double vol = mesh->evalCellVolume(idd);

                    for (std::size_t i = 0; i < nsf; i++){
                        const double *datai = data + scalarIds[i];
                        double *datamappedi = datamapped + scalarIds[i];
                        (*datamappedi) += (*datai) * vol / volmapped;
                    }
                    for (std::size_t i = 0; i < nvf; i++) {
                        for (std::size_t j = 0; j < 3; j++) {
                            const double *datai = data + vectorIds[i][j];
                            double *datamappedi = datamapped + vectorIds[i][j];
                            (*datamappedi) += (*datai) * vol / volmapped;
                        }
                    }
                }
            }
            else if (m_mapper[id].type == mapping::Type::TYPE_REFINEMENT){
                const double *data = fields.data(m_mapper[id].mapped[0]);
                for (std::size_t i = 0; i < nsf; i++){
                    const double *datai = data + scalarIds[i];
                    double *datamappedi = datamapped + scalarIds[i];
                    (*datamappedi) = (*datai);
                }
                for (std::size_t i = 0; i < nvf; i++) {
                    for (std::size_t j = 0; j < 3; j++) {
                        const double *datai = data + vectorIds[i][j];
                        double *datamappedi = datamapped + vectorIds[i][j];
                        (*datamappedi) = (*datai);
                    }
                }
            }
        }//end targetcells

#if BITPIT_ENABLE_MPI
    }
    else{

        //Mapping between partitioned meshes

        //Communicate data and volumes
        std::map<int, std::map<long, std::vector<double> > > datarec;
        std::map<int, std::map<long, double> > volrec;

        communicateField(fields, mesh, datarec, volrec);

        for (long id : *targetCells){
            double *datamapped = mappedFields.data(id);
            if (m_mapper[id].type == mapping::Type::TYPE_RENUMBERING){
                const double *data;
                int rank = m_mapper[id].rank[0];
                if (rank != m_rank){
                    data = datarec[rank][m_mapper[id].mapped[0]].data();
                }
                else{
                    data = fields.data(m_mapper[id].mapped[0]);
                }//endif different rank

                for (std::size_t i = 0; i < nsf; i++){
                    const double *datai = data + scalarIds[i];
                    double *datamappedi = datamapped + scalarIds[i];
                    (*datamappedi) = (*datai);
                }
                for (std::size_t i = 0; i < nvf; i++) {
                    for (std::size_t j = 0; j < 3; j++) {
                        const double *datai = data + vectorIds[i][j];
                        double *datamappedi = datamapped + vectorIds[i][j];
                        (*datamappedi) = (*datai);
                    }
                }
            }
            else if (m_mapper[id].type == mapping::Type::TYPE_COARSENING){
                for (std::size_t i = 0; i < nsf; i++){
                    double *datamappedi = datamapped + scalarIds[i];
                    (*datamappedi) = 0.0;
                }
                for (std::size_t i = 0; i < nvf; i++){
                    for (std::size_t j = 0; j < 3; j++) {
                        double *datamappedi = datamapped + vectorIds[i][j];
                        (*datamappedi) = 0.0;
                    }
                }
                double volmapped = getCellVolume(id);
                int ii=0;
                for (long idd : m_mapper[id].mapped){
                    int rank = m_mapper[id].rank[ii];
                    const double *data = fields.data(idd);
                    double vol = mesh->evalCellVolume(idd);
                    if (rank != m_rank){
                        data = datarec[rank][idd].data();
                        vol = volrec[rank][idd];
                    }
                    else{
                        data = fields.data(idd);
                        vol = mesh->evalCellVolume(idd);
                    }//endif different rank

                    for (std::size_t i = 0; i < nsf; i++){
                        const double *datai = data + scalarIds[i];
                        double *datamappedi = datamapped + scalarIds[i];
                        (*datamappedi) += (*datai) * vol / volmapped;
                    }
                    for (std::size_t i = 0; i < nvf; i++) {
                        for (std::size_t j = 0; j < 3; j++) {
                            const double *datai = data + vectorIds[i][j];
                            double *datamappedi = datamapped + vectorIds[i][j];
                            (*datamappedi) += (*datai) * vol / volmapped;
                        }
                    }
                    ii++;
                }
            }
            else if (m_mapper[id].type == mapping::Type::TYPE_REFINEMENT){
                const double* data;
                int rank = m_mapper[id].rank[0];
                if (rank != m_rank){
                    data = datarec[rank][m_mapper[id].mapped[0]].data();
                }
                else{
                    data = fields.data(m_mapper[id].mapped[0]);
                }//endif different rank
                for (std::size_t i = 0; i < nsf; i++){
                    const double *datai = data + scalarIds[i];
                    double *datamappedi = datamapped + scalarIds[i];
                    (*datamappedi) = (*datai);
                }
                for (std::size_t i = 0; i < nvf; i++) {
                    for (std::size_t j = 0; j < 3; j++) {
                        const double *datai = data + vectorIds[i][j];
                        double *datamappedi = datamapped + vectorIds[i][j];
                        (*datamappedi) = (*datai);
                    }
                }
            }
        }//end targetcells

    }//endif partitioned mesh
#endif

    return mappedFields;

}

/**
 * Map a storage of fields FROM the POD mesh.
 *
 * \param[out] field Output Fields storage
 * \param[in] mesh Pointer to input field mesh
 * \param[in] targetCells Pointer to set of target cells of input field mesh mesh (cells with mapped values as result)
 * \param[in]mappedFields  Fields storage on POD mesh to be mapped.
 * \param[in] scalarIds Ids of scalar fields in PiercedStorage.
 * \param[in] vectorIds Ids of vector fields in PiercedStorage.
 */
void PODVolOctree::mapFieldsFromPOD(PiercedStorage<double> & fields, const VolumeKernel * mesh,
        const std::unordered_set<long> * targetCells, const PiercedStorage<double> & mappedFields,
        const std::vector<std::size_t> &scalarIds, const std::vector<std::array<std::size_t, 3>> &vectorIds)
{

    //Check target cells
    std::unordered_set<long> targetCellsStorage;
    if (!targetCells) {
        for (const Cell &cell : mesh->getCells())
            targetCellsStorage.insert(cell.getId());

        targetCells = &targetCellsStorage;
    }

    // Map data fields from pod mesh
    std::size_t nsf = scalarIds.size();
    std::size_t nvf = vectorIds.size();

    const PiercedStorage<mapping::Info> & m_invmapper = getMapper()->getInverseMapping();

#if BITPIT_ENABLE_MPI
    if ( getMapper()->checkPartition()){
#endif

    for (long id : *targetCells){
        double *data = fields.data(id);
        if (m_invmapper[id].type == mapping::Type::TYPE_RENUMBERING){
            const double *datamapped = mappedFields.data(m_invmapper[id].mapped[0]);
            for (std::size_t i = 0; i < nsf; i++){
                double *datai = data + scalarIds[i];
                const double *datamappedi = datamapped + scalarIds[i];
                (*datai) = (*datamappedi);
            }
            for (std::size_t i = 0; i < nvf; i++) {
                for (std::size_t j = 0; j < 3; j++) {
                    double *datai = data + vectorIds[i][j];
                    const double *datamappedi = datamapped + vectorIds[i][j];
                    (*datai) = (*datamappedi);
                }
            }
        }
        else if (m_invmapper[id].type == mapping::Type::TYPE_COARSENING){
            for (std::size_t i = 0; i < nsf; i++){
                double *datai = data + scalarIds[i];
                (*datai) = 0.0;
            }
            for (std::size_t i = 0; i < nvf; i++){
                for (std::size_t j = 0; j < 3; j++) {
                    double *datai = data + vectorIds[i][j];
                    (*datai) = 0.0;
                }
            }
            double vol = mesh->evalCellVolume(id);
            for (long idd : m_invmapper[id].mapped){
                std::size_t rawIndexIdd = m_meshPOD->getCells().getRawIndex(idd);
                const double * datamapped = mappedFields.data(idd);
                double volmapped = getRawCellVolume(rawIndexIdd);

                for (std::size_t i = 0; i < nsf; i++){
                    double *datai = data + scalarIds[i];
                    const double *datamappedi = datamapped + scalarIds[i];
                    (*datai) += (*datamappedi) * volmapped / vol;
                }
                for (std::size_t i = 0; i < nvf; i++) {
                    for (std::size_t j = 0; j < 3; j++) {
                        double *datai = data + vectorIds[i][j];
                        const double *datamappedi = datamapped + vectorIds[i][j];
                        (*datai) += (*datamappedi) * volmapped / vol;
                    }
                }
            }
        }
        else if (m_invmapper[id].type == mapping::Type::TYPE_REFINEMENT){
            const double *datamapped = mappedFields.data(m_invmapper[id].mapped[0]);
            for (std::size_t i = 0; i < nsf; i++){
                double *datai = data + scalarIds[i];
                const double *datamappedi = datamapped + scalarIds[i];
                (*datai) = (*datamappedi);
            }
            for (std::size_t i = 0; i < nvf; i++) {
                for (std::size_t j = 0; j < 3; j++) {
                    double *datai = data + vectorIds[i][j];
                    const double *datamappedi = datamapped + vectorIds[i][j];
                    (*datai) = (*datamappedi);
                }
            }
        }
    }//end target cells

#if BITPIT_ENABLE_MPI
    }
    else{

        //Mapping between partitioned meshes

        //Communicate data and volumes
        std::map<int, std::map<long, std::vector<double> > > datarec;
        std::map<int, std::map<long, double> > volrec;

        communicateFieldFromPOD(mappedFields, mesh, datarec, volrec);

        for (long id : *targetCells){
            double *data = fields.data(id);
            if (m_invmapper[id].type == mapping::Type::TYPE_RENUMBERING){
                const double *datamapped;
                int rank = m_invmapper[id].rank[0];
                if (rank != m_rank){
                    datamapped = datarec[rank][m_invmapper[id].mapped[0]].data();
                }
                else{
                    datamapped = mappedFields.data(m_invmapper[id].mapped[0]);
                }

                for (std::size_t i = 0; i < nsf; i++){
                    double *datai = data + scalarIds[i];
                    const double *datamappedi = datamapped + scalarIds[i];
                    (*datai) = (*datamappedi);
                }
                for (std::size_t i = 0; i < nvf; i++) {
                    for (std::size_t j = 0; j < 3; j++) {
                        double *datai = data + vectorIds[i][j];
                        const double *datamappedi = datamapped + vectorIds[i][j];
                        (*datai) = (*datamappedi);
                    }
                }
            }
            else if (m_invmapper[id].type == mapping::Type::TYPE_COARSENING){
                for (std::size_t i = 0; i < nsf; i++){
                    double *datai = data + scalarIds[i];
                    (*datai) = 0.0;
                }
                for (std::size_t i = 0; i < nvf; i++){
                    for (std::size_t j = 0; j < 3; j++) {
                        double *datai = data + vectorIds[i][j];
                        (*datai) = 0.0;
                    }
                }
                double vol = mesh->evalCellVolume(id);
                int ii = 0;
                for (long idd : m_invmapper[id].mapped){
                    const double *datamapped;
                    double volmapped;
                    int rank = m_invmapper[id].rank[ii];
                    if (rank != m_rank){
                        datamapped = datarec[rank][idd].data();
                        volmapped = volrec[rank][idd];
                    }
                    else{
                        datamapped = mappedFields.data(idd);
                        volmapped = getCellVolume(idd);
                    }

                    for (std::size_t i = 0; i < nsf; i++){
                        double *datai = data + scalarIds[i];
                        const double *datamappedi = datamapped + scalarIds[i];
                        (*datai) += (*datamappedi) * volmapped / vol;
                    }
                    for (std::size_t i = 0; i < nvf; i++) {
                        for (std::size_t j = 0; j < 3; j++) {
                            double *datai = data + vectorIds[i][j];
                            const double *datamappedi = datamapped + vectorIds[i][j];
                            (*datai) += (*datamappedi) * volmapped / vol;
                        }
                    }
                    ii++;
                }
            }
            else if (m_invmapper[id].type == mapping::Type::TYPE_REFINEMENT){
                const double *datamapped;
                int rank = m_invmapper[id].rank[0];
                if (rank != m_rank){
                    datamapped = datarec[rank][m_invmapper[id].mapped[0]].data();
                }
                else{
                    datamapped = mappedFields.data(m_invmapper[id].mapped[0]);
                }
                for (std::size_t i = 0; i < nsf; i++){
                    double *datai = data + scalarIds[i];
                    const double *datamappedi = datamapped + scalarIds[i];
                    (*datai) = (*datamappedi);
                }
                for (std::size_t i = 0; i < nvf; i++) {
                    for (std::size_t j = 0; j < 3; j++) {
                        double *datai = data + vectorIds[i][j];
                        const double *datamappedi = datamapped + vectorIds[i][j];
                        (*datai) = (*datamappedi);
                    }
                }
            }
        }//end target cells



    }//endif partitioned
#endif


}

/**
 * Map a storage of boolean field TO the POD mesh.
 *
 * \param[in] field Input Fields storage
 * \param[in] mesh Pointer to input field mesh
 * \param[in] targetCells Pointer to set of target cells of POD mesh (cells with mapped values as result)
 * \return Fields storage on target cells of POD mesh.
 */
PiercedStorage<bool> PODVolOctree::mapBoolFieldToPOD(const PiercedStorage<bool> & field, const VolumeKernel * mesh,
        const std::unordered_set<long> * targetCells)
{

    // Map bool field (number of fields = 1) on pod mesh
    PiercedStorage<bool> mappedField(1, &getMesh()->getCells());
    mapBoolFieldToPOD(field, mesh, targetCells, mappedField);
    return mappedField;

}

/**
 * Map a storage of boolean field TO the POD mesh.
 *
 * \param[in] field Input Fields storage
 * \param[in] mesh Pointer to input field mesh
 * \param[in] targetCells Pointer to set of target cells of POD mesh (cells with mapped values as result)
 * \param[out] Fields storage on target cells of POD mesh.
 */
void PODVolOctree::mapBoolFieldToPOD(const PiercedStorage<bool> & field, const VolumeKernel * mesh,
        const std::unordered_set<long> * targetCells, PiercedStorage<bool> & mappedField)
{
    BITPIT_UNUSED(mesh);

    //Check target cells
    std::unordered_set<long> targetCellsStorage;
    if (!targetCells) {
        for (const Cell &cell : getMesh()->getCells())
            targetCellsStorage.insert(cell.getId());

        targetCells = &targetCellsStorage;
    }

    // Map bool field (number of fields = 1) on pod mesh
    const PiercedStorage<mapping::Info> & m_mapper = getMapper()->getMapping();

    mappedField.setStaticKernel(&getMesh()->getCells());
    mappedField.fill(false);

#if BITPIT_ENABLE_MPI
    if ( getMapper()->checkPartition()){
#endif

        for (long id : *targetCells){
            if (m_mapper[id].type == mapping::Type::TYPE_RENUMBERING){
                bool dataB = field.at(m_mapper[id].mapped[0]);
                mappedField.set(id, dataB);
            }
            else if (m_mapper[id].type == mapping::Type::TYPE_COARSENING){
                mappedField.set(id, false);
                bool dataB, dataMappedB = true;
                for (long idd : m_mapper[id].mapped){
                    dataB = field.at(idd);
                    dataMappedB &= dataB;
                }
                mappedField.set(id, dataMappedB);
            }
            else if (m_mapper[id].type == mapping::Type::TYPE_REFINEMENT){
                bool dataB = field.at(m_mapper[id].mapped[0]);
                mappedField.set(id, dataB);
            }
        }

#if BITPIT_ENABLE_MPI
    }
    else{
        //Mapping between partitioned meshes

        //Communicate data
        std::map<int, std::map<long, bool> > dataBrec;

        communicateBoolField(field, dataBrec);

        for (long id : *targetCells){

            if (m_mapper[id].type == mapping::Type::TYPE_RENUMBERING){
                bool dataB;
                int rank = m_mapper[id].rank[0];
                if (rank != m_rank){
                    dataB = dataBrec[rank][m_mapper[id].mapped[0]];
                }
                else{
                    dataB = field.at(m_mapper[id].mapped[0]);
                }//end if other rank
                mappedField.set(id, dataB);
            }
            else if (m_mapper[id].type == mapping::Type::TYPE_COARSENING){
                mappedField.set(id, true);
                bool dataB, dataMappedB = false;
                int ii = 0;
                for (long idd : m_mapper[id].mapped){
                    int rank = m_mapper[id].rank[ii];
                    if (rank != m_rank){
                        dataB = dataBrec[rank][idd];
                    }
                    else{
                        dataB = field.at(idd);
                    }//end if other rank
                    dataMappedB &= dataB;
                    ii++;
                }
                mappedField.set(id, dataMappedB);

            }
            else if (m_mapper[id].type == mapping::Type::TYPE_REFINEMENT){
                bool dataB;
                int rank = m_mapper[id].rank[0];
                if (rank != m_rank){
                    dataB = dataBrec[rank][m_mapper[id].mapped[0]];
                }
                else{
                    dataB = field.at(m_mapper[id].mapped[0]);
                }//end if other rank
                mappedField.set(id, dataB);
            }
        }

    }//end if partitioned
#endif

}

/**
 * Map a set of cells TO the POD mesh by using a pre-computed mapper.
 *
 * \param[in] targetCells Pointer to input cells set
 * return Target cells set of POD mesh.
 */
std::unordered_set<long> PODVolOctree::mapCellsToPOD(const std::unordered_set<long> * targetCells)
{
    std::unordered_set<long> mappedCells;
    const PiercedStorage<mapping::Info> & m_invmapper = getMapper()->getInverseMapping();

#if BITPIT_ENABLE_MPI
    if ( getMapper()->checkPartition()){
#endif

    for (const long & id : *targetCells){
        if (m_invmapper[id].type == mapping::Type::TYPE_RENUMBERING){
            mappedCells.insert(m_invmapper[id].mapped[0]);
        }
        if (m_invmapper[id].type == mapping::Type::TYPE_REFINEMENT){
            mappedCells.insert(m_invmapper[id].mapped[0]);
        }
        if (m_invmapper[id].type == mapping::Type::TYPE_COARSENING){
            for (long idd : m_invmapper[id].mapped)
                mappedCells.insert(idd);
        }
    }//end targetcells

#if BITPIT_ENABLE_MPI
    }
    else{

        //Mapping between partitioned meshes

        //Inset local POD cells and recover POD cells on other partitions by inverse mapping
        std::map<int, std::unordered_set<long> > sendPODcells;

        for (const long & id : *targetCells){
            if (m_invmapper[id].type == mapping::Type::TYPE_RENUMBERING){
                int rank = m_invmapper[id].rank[0];
                if (rank != m_rank){
                    sendPODcells[rank].insert(m_invmapper[id].mapped[0]);
                }
                else{
                    mappedCells.insert(m_invmapper[id].mapped[0]);
                }
            }
            if (m_invmapper[id].type == mapping::Type::TYPE_REFINEMENT){
                int rank = m_invmapper[id].rank[0];
                if (rank != m_rank){
                    sendPODcells[rank].insert(m_invmapper[id].mapped[0]);
                }
                else{
                    mappedCells.insert(m_invmapper[id].mapped[0]);
                }
            }
            if (m_invmapper[id].type == mapping::Type::TYPE_COARSENING){
                int ii=0;
                for (long idd : m_invmapper[id].mapped){
                    int rank = m_invmapper[id].rank[ii];
                    if (rank != m_rank){
                        sendPODcells[rank].insert(idd);
                    }
                    else{
                        mappedCells.insert(idd);
                    }
                    ii++;
                }
            }
        }//end targetcells


        //Communicate cells to send
        //build send buffers
        DataCommunicator dataCommunicator(m_communicator);
        std::size_t bytes = uint8_t(sizeof(long));
        for (std::pair<int, std::unordered_set<long> > val : sendPODcells){
            std::size_t ncells = val.second.size();
            int rank = val.first;
            //set size
            std::size_t buffSize = ncells * bytes;
            dataCommunicator.setSend(rank,buffSize);
            //fill buffer with octants
            SendBuffer &sendBuffer = dataCommunicator.getSendBuffer(rank);
            for (const long & ID : val.second){
                sendBuffer << ID;
            }
        }

        dataCommunicator.discoverRecvs();
        dataCommunicator.startAllRecvs();
        dataCommunicator.startAllSends();

        vector<int> recvRanks = dataCommunicator.getRecvRanks();
        std::sort(recvRanks.begin(),recvRanks.end());

        for (int rank : recvRanks){
            dataCommunicator.waitRecv(rank);
            RecvBuffer & recvBuffer = dataCommunicator.getRecvBuffer(rank);
            long nof = recvBuffer.getSize() / bytes;
            for (long ii = 0; ii < nof; ii++){
                long ID;
                recvBuffer >> ID;
                mappedCells.insert(ID);
            }
        }

        dataCommunicator.waitAllSends();

    }//endif partitioned
#endif

    return mappedCells;

}

void PODVolOctree::adaptMeshToMesh(VolumeKernel* meshToAdapt, VolumeKernel * meshReference)
{
    BITPIT_UNUSED(meshToAdapt);

    computeMapper(meshReference, false);

    bool adapt = true;

    while(adapt){

        adapt = false;

        const PiercedStorage<mapping::Info> & m_mapper = getMapper()->getMapping();

        for (Cell & cell : getMesh()->getCells()){
            long id = cell.getId();
            if (m_mapper[id].type == mapping::Type::TYPE_COARSENING){
                getMesh()->markCellForRefinement(id);
                adapt = true;
            }
        }

#if BITPIT_ENABLE_MPI
        MPI_Allreduce(MPI_IN_PLACE, &adapt, 1, MPI_C_BOOL, MPI_LOR, m_communicator);
#endif

        std::vector<adaption::Info> infoAdapt = getMesh()->adaptionPrepare(true);

        getMapper()->adaptionPrepare(infoAdapt);

        infoAdapt = getMesh()->adaptionAlter(true);

        getMapper()->adaptionAlter(infoAdapt, true, false);

        getMesh()->adaptionCleanup();

        getMapper()->adaptionCleanup();

    }

}

# if BITPIT_ENABLE_MPI

void PODVolOctree::communicatePODField(const pod::PODField & field, std::map<int, std::map<long, bool> > & dataBrec, std::map<int, std::map<long, std::vector<double> > > & dataSrec, std::map<int, std::map<long, std::vector<std::array<double,3> > > > & dataVrec, std::map<int, std::map<long, double> > & volrec)
{

    {
        // Size of fields
        std::size_t nsf = field.scalar->getFieldCount();
        std::size_t nvf = field.vector->getFieldCount();

        std::map<int, std::vector<long> > rankIDrec = m_meshmap->getReceivedMappedID();
        std::map<int, std::vector<long> > rankIDsend = m_meshmap->getSentMappedID();

        //build send buffers
        DataCommunicator dataCommunicator(m_communicator);
        std::size_t bytes = uint8_t(sizeof(bool) + (nsf+nvf+1)*sizeof(double));
        for (std::pair<const int, std::vector<long int> > val : rankIDsend){
            int rank = val.first;
            //set size
            std::size_t buffSize = val.second.size() * bytes;
            dataCommunicator.setSend(rank,buffSize);
            //fill buffer with octants
            SendBuffer &sendBuffer = dataCommunicator.getSendBuffer(rank);
            for (long & ID : val.second){
                sendBuffer << field.mask->at(ID);
                for (std::size_t i=0; i<nsf; i++)
                    sendBuffer << field.scalar->at(ID,i);
                for (std::size_t i=0; i<nvf; i++)
                    sendBuffer << field.vector->at(ID,i);
                sendBuffer << field.mesh->evalCellVolume(ID);
            }
        }

        dataCommunicator.startAllRecvs();
        dataCommunicator.startAllSends();

        for (std::pair<const int, std::vector<long int> > val : rankIDrec){
            int rank = val.first;
            dataCommunicator.waitRecv(rank);
            RecvBuffer & recvBuffer = dataCommunicator.getRecvBuffer(rank);
            long bufferSize = recvBuffer.getSize();
            for (long & ID : val.second){
                recvBuffer >> dataBrec[rank][ID];
                for (std::size_t i=0; i<nsf; i++){
                    double data;
                    recvBuffer >> data;
                    dataSrec[rank][ID].push_back(data);
                }
                for (std::size_t i=0; i<nvf; i++){
                    std::array<double,3> data;
                    recvBuffer >> data[0];
                    recvBuffer >> data[1];
                    recvBuffer >> data[2];
                    dataVrec[rank][ID].push_back(data);
                }
                recvBuffer >> volrec[rank][ID];
            }
        }
        dataCommunicator.waitAllSends();

    }

}

void PODVolOctree::communicatePODFieldFromPOD(const pod::PODField & field, std::map<int, std::map<long, bool> > & dataBrec, std::map<int, std::map<long, std::vector<double> > > & dataSrec, std::map<int, std::map<long, std::vector<std::array<double,3> > > > & dataVrec, std::map<int, std::map<long, double> > & volrec)
{

    {
        // Size of fields
        std::size_t nsf = field.scalar->getFieldCount();
        std::size_t nvf = field.vector->getFieldCount();

        std::map<int, std::vector<long> > rankIDsend = m_meshmap->getSentReferenceID();

        //build send buffers
        DataCommunicator dataCommunicator(m_communicator);
        std::size_t bytes = uint8_t(sizeof(long) + sizeof(bool) + (nsf+nvf+1)*sizeof(double));
        for (std::pair<const int, std::vector<long int> > val : rankIDsend){
            int rank = val.first;
            //set size
            std::size_t buffSize = val.second.size() * bytes;
            dataCommunicator.setSend(rank,buffSize);
            //fill buffer with octants
            SendBuffer &sendBuffer = dataCommunicator.getSendBuffer(rank);
            for (long & ID : val.second){
                sendBuffer << ID;
                sendBuffer << field.mask->at(ID);
                for (std::size_t i=0; i<nsf; i++)
                    sendBuffer << field.scalar->at(ID,i);
                for (std::size_t i=0; i<nvf; i++)
                    sendBuffer << field.vector->at(ID,i);
                sendBuffer << field.mesh->evalCellVolume(ID);
            }
        }

        dataCommunicator.discoverRecvs();
        dataCommunicator.startAllRecvs();
        dataCommunicator.startAllSends();

        vector<int> recvRanks = dataCommunicator.getRecvRanks();
        std::sort(recvRanks.begin(),recvRanks.end());

        for (int rank : recvRanks){
            dataCommunicator.waitRecv(rank);
            RecvBuffer & recvBuffer = dataCommunicator.getRecvBuffer(rank);
            long nof = recvBuffer.getSize() / bytes;
            for (long ii = 0; ii < nof; ii++){
                long ID;
                recvBuffer >> ID;
                recvBuffer >> dataBrec[rank][ID];
                for (std::size_t i=0; i<nsf; i++){
                    double data;
                    recvBuffer >> data;
                    dataSrec[rank][ID].push_back(data);
                }
                for (std::size_t i=0; i<nvf; i++){
                    std::array<double,3> data;
                    recvBuffer >> data[0];
                    recvBuffer >> data[1];
                    recvBuffer >> data[2];
                    dataVrec[rank][ID].push_back(data);
                }
                recvBuffer >> volrec[rank][ID];
            }
        }

        dataCommunicator.waitAllSends();

    }

}

void PODVolOctree::communicateBoolField(const PiercedStorage<bool> & field, std::map<int, std::map<long, bool> > & dataBrec)
{

    std::map<int, std::vector<long> > rankIDrec = m_meshmap->getReceivedMappedID();
    std::map<int, std::vector<long> > rankIDsend = m_meshmap->getSentMappedID();

    //build send buffers
    DataCommunicator dataCommunicator(m_communicator);
    std::size_t bytes = uint8_t(sizeof(bool));
    for (std::pair<const int, std::vector<long int> > val : rankIDsend){
        int rank = val.first;
        //set size
        std::size_t buffSize = val.second.size() * bytes;
        dataCommunicator.setSend(rank,buffSize);
        //fill buffer with octants
        SendBuffer &sendBuffer = dataCommunicator.getSendBuffer(rank);
        for (long & ID : val.second){
            sendBuffer << field.at(ID);
        }
    }

    dataCommunicator.startAllRecvs();
    dataCommunicator.startAllSends();

    for (std::pair<const int, std::vector<long int> > val : rankIDrec){
        int rank = val.first;
        dataCommunicator.waitRecv(rank);
        RecvBuffer & recvBuffer = dataCommunicator.getRecvBuffer(rank);
        long bufferSize = recvBuffer.getSize();
        for (long & ID : val.second){
            recvBuffer >> dataBrec[rank][ID];
        }
    }

    dataCommunicator.waitAllSends();

}

void PODVolOctree::communicateField(const PiercedStorage<double> & field, const VolumeKernel * mesh, std::map<int, std::map<long, std::vector<double> > > & datarec, std::map<int, std::map<long, double> > & volrec)
{

    std::size_t nf = field.getFieldCount();

    std::map<int, std::vector<long> > rankIDrec = m_meshmap->getReceivedMappedID();
    std::map<int, std::vector<long> > rankIDsend = m_meshmap->getSentMappedID();

    //build send buffers
    DataCommunicator dataCommunicator(m_communicator);
    std::size_t bytes = uint8_t((nf+1)*sizeof(double));
    for (std::pair<const int, std::vector<long int> > val : rankIDsend){
        int rank = val.first;
        //set size
        std::size_t buffSize = val.second.size() * bytes;
        dataCommunicator.setSend(rank,buffSize);
        //fill buffer with octants
        SendBuffer &sendBuffer = dataCommunicator.getSendBuffer(rank);
        for (long & ID : val.second){
            for (std::size_t ifield=0; ifield<nf; ifield++)
                sendBuffer << field.at(ID, ifield);
            sendBuffer << mesh->evalCellVolume(ID);
        }
    }

    dataCommunicator.startAllRecvs();
    dataCommunicator.startAllSends();

    for (std::pair<const int, std::vector<long int> > val : rankIDrec){
        int rank = val.first;
        dataCommunicator.waitRecv(rank);
        RecvBuffer & recvBuffer = dataCommunicator.getRecvBuffer(rank);
        long bufferSize = recvBuffer.getSize();
        for (long & ID : val.second){
            datarec[rank][ID].resize(nf);
            for (std::size_t ifield=0; ifield<nf; ifield++)
                recvBuffer >> datarec[rank][ID][ifield];
            recvBuffer >> volrec[rank][ID];
        }
    }

    dataCommunicator.waitAllSends();

}

void PODVolOctree::communicateFieldFromPOD(const PiercedStorage<double> & field, const VolumeKernel * mesh, std::map<int, std::map<long, std::vector<double> > > & datarec, std::map<int, std::map<long, double> > & volrec)
{

    std::size_t nf = field.getFieldCount();

    std::map<int, std::vector<long> > rankIDsend = m_meshmap->getSentReferenceID();

    //build send buffers
    DataCommunicator dataCommunicator(m_communicator);
    std::size_t bytes = uint8_t((nf+1)*sizeof(double)+sizeof(long));
    for (std::pair<const int, std::vector<long int> > val : rankIDsend){
        int rank = val.first;
        //set size
        std::size_t buffSize = val.second.size() * bytes;
        dataCommunicator.setSend(rank,buffSize);
        //fill buffer with octants
        SendBuffer &sendBuffer = dataCommunicator.getSendBuffer(rank);
        for (long & ID : val.second){
            sendBuffer << ID;
            for (std::size_t ifield=0; ifield<nf; ifield++)
                sendBuffer << field.at(ID, ifield);
            sendBuffer << mesh->evalCellVolume(ID);
        }
    }

    dataCommunicator.discoverRecvs();
    dataCommunicator.startAllRecvs();
    dataCommunicator.startAllSends();

    vector<int> recvRanks = dataCommunicator.getRecvRanks();
    std::sort(recvRanks.begin(),recvRanks.end());

    for (int rank : recvRanks){
        dataCommunicator.waitRecv(rank);
        RecvBuffer & recvBuffer = dataCommunicator.getRecvBuffer(rank);
        long nof = recvBuffer.getSize() / bytes;
        for (long ii = 0; ii < nof; ii++){
            long ID;
            recvBuffer >> ID;
            datarec[rank][ID].resize(nf);
            for (std::size_t ifield=0; ifield<nf; ifield++)
                recvBuffer >> datarec[rank][ID][ifield];
            recvBuffer >> volrec[rank][ID];
        }
    }

    dataCommunicator.waitAllSends();

}

#endif


}
