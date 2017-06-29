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

# include <vector>
# include <stack>

# if BITPIT_ENABLE_MPI
# include <mpi.h>
# include "bitpit_communications.hpp"
# endif

# include "bitpit_operators.hpp"
# include "bitpit_CG.hpp"
# include "bitpit_patchkernel.hpp"

# include "levelSetKernel.hpp"
# include "levelSetObject.hpp"

namespace bitpit {

/*!
	@interface LevelSetObject
	@ingroup levelset
	@brief Interface class for all objects with respect to whom the levelset function may be computed.
*/

/*!
 * Destructor
 */
LevelSetObject::~LevelSetObject( ){
}

/*!
 * Constructor
 * @param[in] id id assigned to object
 */
LevelSetObject::LevelSetObject(int id) : m_id(id){
    m_kernelPtr=nullptr;
}

/*!
 * Sets the kernel for the object
 * @param[in] kernel is the LevelSetKernel
 */
void LevelSetObject::setKernel(LevelSetKernel *kernel) {
    m_kernelPtr = kernel;
}

/*!
 * Get the id 
 * @return id of the object
 */
int LevelSetObject::getId( ) const {
    return m_id ;
}

/*!
 * If the levelset is primary (e.g. of a surface triangulation) or not (e.g. derived by boolean operations between two levelsets)
 * @return if object is primary
 */
bool LevelSetObject::isPrimary( ) const {
    return true;
}

/*!
 * Get the object and part id of projection point
 * @param[in] i cell index
 * @return pair containing object and part id 
 */
int LevelSetObject::getPart(const long &i) const {
    BITPIT_UNUSED(i) ;
    return levelSetDefaults::PART ;
}

/*!
 * Get the sign of the levelset function
 * @param[in] i cell index
 * @return sign of levelset
 */
short LevelSetObject::getSign(const long &i)const{
    return ( static_cast<short>(sign(getLS(i) )) );
}

/*!
 * Propgates the sign to levelset function throughout the grid
 */
void LevelSetObject::propagateSign(){
}

/*!
 * If cell centroid lies within the narrow band and hence levelset is computet exactly
 * @param[in] i cell index
 * @return true/false if the centroid is in narrow band
 */
bool LevelSetObject::isInNarrowBand(const long &i)const{
    return ( std::abs(getLS(i)) <= m_RSearch );
}

/*!
 * Get the current size of the narrow band.
 * @return size of the current narrow band
 */
double LevelSetObject::getSizeNarrowBand()const{
    return m_RSearch;
}

/*!
 * Manually set the size of the narrow band.
 * @param[in] r size of the narrow band.
 */
void LevelSetObject::setSizeNarrowBand(double r){
    m_RSearch = r;
}

/*!
 * Check if cell intersects the surface
 * @param[in] i cell index
 * @return indicator regarding intersection
 */
LevelSetIntersectionStatus LevelSetObject::intersectSurface(const long &i) const{
    if( !isInNarrowBand(i)){
        return LevelSetIntersectionStatus::FALSE;
    }

    double circumcircle = m_kernelPtr->computeCellCircumcircle(i) ;
    if(std::abs(getLS(i)) > circumcircle){
        return LevelSetIntersectionStatus::FALSE;
    }

    double incircle = m_kernelPtr->computeCellIncircle(i) ;
    if(std::abs(getLS(i)) <= incircle){
        return LevelSetIntersectionStatus::TRUE;
    }

    return LevelSetIntersectionStatus::CLOSE;
}

/*!
 * Returns the characterstic size at the support
 * @param[in] i cell index
 * @return feature size
 */
double LevelSetObject::getSurfaceFeatureSize(const long &i) const{
    BITPIT_UNUSED(i);
    return (- levelSetDefaults::SIZE);
}

/*!
 * Returns the minimum surface feature size
 * @return feature size
 */
double LevelSetObject::getMinSurfaceFeatureSize() const{
    return (- levelSetDefaults::SIZE);
}

/*!
 * Returns the maximum surface feature size
 * @return feature size
 */
double LevelSetObject::getMaxSurfaceFeatureSize() const{
    return (- levelSetDefaults::SIZE);
}

/*!
 * Calculates the value and gradient of the levelset function within the narrow band
 * @param[in] RSeach size of the narrow band.
 * @param[in] signd if signed distances should be calculted
 */
void LevelSetObject::computeLSInNarrowBand(const double &RSearch, const bool &signd){
    BITPIT_UNUSED(RSearch);
    BITPIT_UNUSED(signd);
}

/*!
 * Updates the size of the narrow band after mesh adaption
 * @param[in] mapper information regarding mesh adaption
 * @return size of narrow band
 */
double LevelSetObject::updateSizeNarrowBand(const std::vector<adaption::Info> &mapper){
    BITPIT_UNUSED(mapper);
    return getSizeNarrowBand();
}

/*!
 * Updates the value and gradient of the levelset function within the narrow band
 * @param[in] RSeach size of the narrow band.
 * @param[in] signd if signed distances should be calculted
 */
void LevelSetObject::updateLSInNarrowBand(const std::vector<adaption::Info> &mapper, const double &RSearch, const bool &signd){
    BITPIT_UNUSED(mapper);
    BITPIT_UNUSED(RSearch);
    BITPIT_UNUSED(signd);
}

/*! 
 * Deletes non-existing items and items outside the narrow band after grid adaption.
 * @param[in] mapper mapping info
 */
void LevelSetObject::clearAfterMeshAdaption( const std::vector<adaption::Info> &mapper ){
    _clearAfterMeshAdaption( mapper ) ;
}

/*!
 * Clears data structure after mesh modification
 * @param[in] mapper mapper describing mesh modifications
 */
void LevelSetObject::_clearAfterMeshAdaption( const std::vector<adaption::Info> &mapper ){
    BITPIT_UNUSED(mapper) ;
}

/*! 
 * Deletes items outside the narrow band after grid adaption.
 * @param[in] newRSearch new size of narrow band
 */
void LevelSetObject::filterOutsideNarrowBand( double newRSearch ){
    _filterOutsideNarrowBand(newRSearch) ;
}

/*!
 * Clears data structure outside narrow band
 * @param[in] search size of narrow band
 */
void LevelSetObject::_filterOutsideNarrowBand( double search ){
    BITPIT_UNUSED(search) ;
}

/*! 
 * Clears all levelset information
 */
void LevelSetObject::clear( ){
    _clear() ;
}

/*! 
 * Clears all levelset information stored in derived class
 */
void LevelSetObject::_clear( ){
}

/*!
 * Writes LevelSetObject to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetObject::dump( std::ostream &stream ){
    utils::binary::write(stream, m_id) ;
    utils::binary::write(stream, m_RSearch);
    _dump(stream) ;
}

/*!
 * Writes LevelSetObject to stream in binary format
 * @param[in] stream output stream
 */
void LevelSetObject::_dump( std::ostream &stream ){
    BITPIT_UNUSED(stream);
}

/*!
 * Reads LevelSetObject from stream in binary format
 * @param[in] stream output stream
 */
void LevelSetObject::restore( std::istream &stream ){
    utils::binary::read(stream, m_id) ;
    utils::binary::read(stream, m_RSearch);
    _restore(stream) ;
}

/*!
 * Reads LevelSetObject from stream in binary format
 * @param[in] stream output stream
 */
void LevelSetObject::_restore( std::istream &stream ){
    BITPIT_UNUSED(stream);
}

#if BITPIT_ENABLE_MPI

/*!
 * Checks if MPI communicator is available in underlying mesh.
 * If available, MPI communicator is retreived from mesh (and duplicated if necessary) and parallel processing can be done.
 * If not serial processing is necessary
 * @return true if parallel
 */
bool LevelSetObject::assureMPI(){

    return(m_kernelPtr->assureMPI() ) ;
}

/*!
 * Update of ghost cell;
 */
void LevelSetObject::exchangeGhosts(){

    std::unordered_map<int,std::vector<long>> &sendList =  m_kernelPtr->getMesh()->getGhostExchangeSources() ;
    std::unordered_map<int,std::vector<long>> &recvList =  m_kernelPtr->getMesh()->getGhostExchangeTargets() ;

    communicate(sendList,recvList);
}

/*!
 * communicates data structures of kernel and objects.
 * If mapper!=NULL, clearAfterMeshAdaption of kernel and objects will be called between send and receive.
 * @param[in] sendList list of elements to be send
 * @param[in] recvList list of elements to be received
 * @param[in] mapper mapper containing mesh modifications
 */
void LevelSetObject::communicate( std::unordered_map<int,std::vector<long>> &sendList, std::unordered_map<int,std::vector<long>> &recvList, std::vector<adaption::Info> const *mapper){

    if( assureMPI() ){

        MPI_Comm meshComm = m_kernelPtr->getCommunicator() ;

        DataCommunicator    sizeCommunicator(meshComm) ; 
        DataCommunicator    dataCommunicator(meshComm) ; 

        int offset = 2*getId() ;
        sizeCommunicator.setTag(offset+0) ;
        dataCommunicator.setTag(offset+1) ;

        /* Start receiving of first communications which contain the sizes
         * of the communications of the exchanged data
        */
	    for (const auto entry : sendList) {
            int rank = entry.first;

            sizeCommunicator.setRecv( rank, sizeof(long) ) ;
            sizeCommunicator.startRecv(rank);
        }

        /* Fill the databuffer with its content from the LevelSetObject 
         * base class and specific derived class. Fill the sizebuffer
         * with the size of the databuffer. Start sending both buffers
         */
	    for (const auto entry : sendList) {

	    	// Get the send buffer
	    	int rank = entry.first;

            sizeCommunicator.setSend(rank, sizeof(long) ) ;
            dataCommunicator.setSend(rank, 0 ) ;

            SendBuffer &sizeBuffer = sizeCommunicator.getSendBuffer(rank);
            SendBuffer &dataBuffer = dataCommunicator.getSendBuffer(rank);

            // Write data in databuffer
            writeCommunicationBuffer( entry.second, dataBuffer ) ;

            // Write size of databuffer in sizebuffer
            sizeBuffer << dataBuffer.capacity() ;

            // Start sending
            sizeCommunicator.startSend(rank);
            dataCommunicator.startSend(rank);
	    }

        
        /* Check which size communications have arrived. For those which are
         * available read the size of the databuffer and start receiving.
         * Proceed once all sizeBuffers have arrived.
         */
        int nCompletedRecvs(0);
        long dataSize;

        while (nCompletedRecvs < sizeCommunicator.getRecvCount()) {
            int rank = sizeCommunicator.waitAnyRecv();
            RecvBuffer &sizeBuffer = sizeCommunicator.getRecvBuffer(rank);

            sizeBuffer >> dataSize ;

            dataCommunicator.setRecv(rank,dataSize) ;
            dataCommunicator.startRecv(rank) ;

            ++nCompletedRecvs;
        }
        

        /* In case of mesh adaption, the ids of the moved items must be freed 
         * in order to avoid conflicts when ids are recycled
         */
        if( mapper!=NULL){
            clearAfterMeshAdaption(*mapper);
        }

        /* Check which data communications have arrived. For those which are
         * available start reading the databuffer into the data structure of
         * LevelSetObject and its derived classes.
         */
        nCompletedRecvs = 0;
        while (nCompletedRecvs < dataCommunicator.getRecvCount()) {
            int rank = dataCommunicator.waitAnyRecv();

            RecvBuffer &dataBuffer = dataCommunicator.getRecvBuffer(rank);
            readCommunicationBuffer( recvList[rank], dataBuffer ) ;

            ++nCompletedRecvs;
        }
        

        sizeCommunicator.waitAllSends();
        sizeCommunicator.finalize() ;

        dataCommunicator.waitAllSends();
        dataCommunicator.finalize() ;

    }
}

/*!
 * Flushing of data to communication buffers for partitioning
 * @param[in] sendList list of cells to be sent
 * @param[in,out] dataBuffer buffer for second communication containing data
 */
void LevelSetObject::writeCommunicationBuffer( const std::vector<long> &sendList, SendBuffer &dataBuffer ){
    _writeCommunicationBuffer( sendList, dataBuffer) ;
    dataBuffer.squeeze( ) ;
}

/*!
 * Flushing of data to communication buffers for partitioning
 * @param[in] sendList list of cells to be sent
 * @param[in,out] dataBuffer buffer for second communication containing data
 */
void LevelSetObject::_writeCommunicationBuffer( const std::vector<long> &sendList, SendBuffer &dataBuffer ){
    BITPIT_UNUSED(sendList) ;
    BITPIT_UNUSED(dataBuffer) ;
}

/*!
 * Processing of communication buffer into data structure
 * @param[in] recvList list of cells to be received
 * @param[in,out] dataBuffer buffer containing the data
 */
void LevelSetObject::readCommunicationBuffer( const std::vector<long> &recvList, RecvBuffer &dataBuffer ){
    _readCommunicationBuffer(recvList, dataBuffer) ;
}

/*!
 * Processing of communication buffer into data structure
 * @param[in] recvList list of cells to be received
 * @param[in,out] dataBuffer buffer containing the data
 */
void LevelSetObject::_readCommunicationBuffer( const std::vector<long> &recvList, RecvBuffer &dataBuffer ){
    BITPIT_UNUSED(recvList) ;
    BITPIT_UNUSED(dataBuffer) ;
}

#endif 

}
