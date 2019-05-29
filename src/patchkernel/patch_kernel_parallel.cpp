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
#if BITPIT_ENABLE_MPI==1

// ========================================================================== //
// INCLUDES                                                                   //
// ========================================================================== //
#include <mpi.h>
#include <chrono>
#include <unordered_set>

#include "bitpit_SA.hpp"

#include "patch_kernel.hpp"

// ========================================================================== //
// NAMESPACES                                                                 //
// ========================================================================== //
using namespace std;
using namespace chrono;

namespace bitpit {

/*!
	Sets the MPI communicator to be used for parallel communications.

	\param communicator is the communicator to be used for parallel
	communications.
*/
void PatchKernel::setCommunicator(MPI_Comm communicator)
{
	// Communication can be set just once
	if (isCommunicatorSet()) {
		throw std::runtime_error ("Patch communicator can be set just once");
	}

	// The communicator has to be valid
	if (communicator == MPI_COMM_NULL) {
		throw std::runtime_error ("Patch communicator is not valid");
	}

	// Creat a copy of the user-specified communicator
	//
	// No library routine should use MPI_COMM_WORLD as the communicator;
	// instead, a duplicate of a user-specified communicator should always
	// be used.
	MPI_Comm_dup(communicator, &m_communicator);

	// Get MPI information
	MPI_Comm_size(m_communicator, &m_nProcessors);
	MPI_Comm_rank(m_communicator, &m_rank);

	// Set parallel data for the VTK output
	if (m_nProcessors > 1) {
		m_vtk.setParallel(m_nProcessors, m_rank);
	}
}

/*!
	Checks if the communicator to be used for parallel communications has
	already been set.

	\result Returns true if the communicator has been set, false otherwise.
*/
bool PatchKernel::isCommunicatorSet() const
{
	return (getCommunicator() != MPI_COMM_NULL);
}

/*!
	Gets the MPI communicator associated to the patch

	\return The MPI communicator associated to the patch.
*/
const MPI_Comm & PatchKernel::getCommunicator() const
{
	return m_communicator;
}

/*!
	Frees the MPI communicator associated to the patch
*/
void PatchKernel::freeCommunicator()
{
	if (!isCommunicatorSet()) {
		return;
	}

	int finalizedCalled;
	MPI_Finalized(&finalizedCalled);
	if (finalizedCalled) {
		return;
	}

	MPI_Comm_free(&m_communicator);
}

/*!
	Gets the MPI rank associated to the patch

	\return The MPI rank associated to the patch.
*/
int PatchKernel::getRank() const
{
	return m_rank;
}

/*!
	Gets the MPI processors in the communicator associated to the patch

	\return The MPI processors in the communicator associated to the patch
*/
int PatchKernel::getProcessorCount() const
{
	return m_nProcessors;
}

/*!
	Sets the size, expressed in number of layers, of the ghost cells halo.

	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
*/
void PatchKernel::setHaloSize(std::size_t haloSize)
{
	if (isPartitioned()) {
		throw std::runtime_error ("Halo size can only be set before partitionig the patch.");
	}

	std::size_t maxHaloSize = _getMaxHaloSize();
	if (haloSize > maxHaloSize) {
		throw std::runtime_error ("Halo size exceeds the maximum allowed value.");
	}

	if (m_haloSize == haloSize) {
        return;
    }

	m_haloSize = haloSize;

	_setHaloSize(haloSize);

	if (isPartitioned()) {
		updateGhostExchangeInfo();
	}
}

/*!
	Gets the size, expressed in number of layers, of the ghost cells halo.

	\result The size, expressed in number of layers, of the ghost cells halo.
*/
std::size_t PatchKernel::getHaloSize() const
{
	return m_haloSize;
}

/*!
	Gets the maximum allowed size, expressed in number of layers, of the ghost
	cells halo.

	\result The maximum allowed size, expressed in number of layers, of the
	ghost cells halo.
*/
std::size_t PatchKernel::_getMaxHaloSize()
{
	return 1;
}

/*!
	Internal function to set the size, expressed in number of layers, of the
	ghost cells halo.

	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
*/
void PatchKernel::_setHaloSize(std::size_t haloSize)
{
	BITPIT_UNUSED(haloSize);
}

/*!
	Converts an internal cell to a ghost cell.

	\param[in] id is the index of the cell
	\param[in] ownerRank is the owner of the cell
*/
PatchKernel::CellIterator PatchKernel::moveInternal2Ghost(const long &id, int ownerRank)
{
	if (!isExpert()) {
		return m_cells.end();
	}

	// Swap the element with the last internal cell
	if (id != m_lastInternalId) {
		m_cells.swap(id, m_lastInternalId);
	}

	// Get the iterator pointing to the updated position of the element
	CellIterator iterator = m_cells.find(id);

	// Update the interior flag
	iterator->setInterior(false);

	// Update cell counters
	--m_nInternals;
	++m_nGhosts;

	// Update the last internal and first ghost markers
	m_firstGhostId = id;
	if (m_nInternals == 0) {
		m_lastInternalId = Cell::NULL_ID;
	} else {
		m_lastInternalId = m_cells.getSizeMarker(m_nInternals - 1, Cell::NULL_ID);
	}

	// Set ghost owner
	setGhostOwner(id, ownerRank);

	// Return the iterator to the new position
	return iterator;
}

/*!
	Converts a ghost cell to an internal cell.

	\param[in] id is the index of the cell
*/
PatchKernel::CellIterator PatchKernel::moveGhost2Internal(const long &id)
{
	if (!isExpert()) {
		return m_cells.end();
	}

	// Swap the cell with the first ghost
	if (id != m_firstGhostId) {
		m_cells.swap(id, m_firstGhostId);
	}

	// Get the iterator pointing to the updated position of the element
	CellIterator iterator = m_cells.find(id);

	// Update the interior flag
	iterator->setInterior(true);

	// Update cell counters
	++m_nInternals;
	--m_nGhosts;

	// Update the last internal and first ghost markers
	m_lastInternalId = id;
	if (m_nGhosts == 0) {
		m_firstGhostId = Cell::NULL_ID;
	} else {
		CellIterator firstGhostIterator = iterator;
		++firstGhostIterator;
		m_firstGhostId = firstGhostIterator->getId();
	}

	// Unset ghost owner
	unsetGhostOwner(id);

	// Return the iterator to the new position
	return iterator;
}

/*!
	Gets the number of ghost cells in the patch.

	\return The number of ghost cells in the patch
*/
long PatchKernel::getGhostCount() const
{
	return m_nGhosts;
}

/*!
	Gets a reference to the first ghost cell.

	\return A reference to the first ghost cell.
*/
Cell & PatchKernel::getFirstGhost()
{
	return m_cells[m_firstGhostId];
}

/*!
	Gets a constant reference to the first ghost cell.

	\return A constant reference to the first ghost cell.
*/
const Cell & PatchKernel::getFirstGhost() const
{
	return m_cells[m_firstGhostId];
}

/*!
	Adds the specified cell to the patch.

	\param source is the cell that will be added
	\param rank is the rank that owns the cell that will be added
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(const Cell &source, int rank, long id)
{
	if (id < 0) {
		id = generateCellId();
	}

	Cell cell = source;

	return addCell(std::move(cell), rank, id);
}

/*!
	Adds the specified cell to the patch.

	\param source is the cell that will be added
	\param rank is the rank that owns the cell that will be added
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, the id of the source will be used
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(Cell &&source, int rank, long id)
{
	if (id < 0) {
		id = source.getId();
	}

	int connectSize = source.getConnectSize();
	std::unique_ptr<long[]> connectStorage = std::unique_ptr<long[]>(new long[connectSize]);
	if (!source.hasInfo()){
		std::copy(source.getConnect(), source.getConnect() + connectSize, connectStorage.get());
	}

	CellIterator iterator = addCell(source.getType(), std::move(connectStorage), rank, id);

	Cell &cell = (*iterator);
	id = cell.getId();
	cell = std::move(source);
	cell.setId(id);

	return iterator;
}

/*!
	Adds a new cell with the specified id and type.

	\param type is the type of the cell
	\param rank is the rank that owns the cell that will be added
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementType type, int rank, long id)
{
	std::unique_ptr<long[]> connectStorage;
	if (ReferenceElementInfo::hasInfo(type)) {
		int connectSize = ReferenceElementInfo::getInfo(type).nVertices;
		connectStorage = std::unique_ptr<long[]>(new long[connectSize]);
	} else {
		connectStorage = std::unique_ptr<long[]>(nullptr);
	}

	return addCell(type, std::move(connectStorage), rank, id);
}

/*!
	Adds a new cell with the specified id, type, and connectivity.

	\param type is the type of the cell
	\param connectivity is the connectivity of the cell
	\param rank is the rank that owns the cell that will be added
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementType type, const std::vector<long> &connectivity,
											   int rank, long id)
{
	int connectSize = connectivity.size();
	std::unique_ptr<long[]> connectStorage = std::unique_ptr<long[]>(new long[connectSize]);
	std::copy(connectivity.data(), connectivity.data() + connectSize, connectStorage.get());

	return addCell(type, std::move(connectStorage), rank, id);
}

/*!
	Adds a new cell with the specified id, type, and connectivity.

	\param type is the type of the cell
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
	\param rank is the rank that owns the cell that will be added
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the added cell.
*/
PatchKernel::CellIterator PatchKernel::addCell(ElementType type, std::unique_ptr<long[]> &&connectStorage,
											   int rank, long id)
{
	if (!isExpert()) {
		return cellEnd();
	}

	if (id < 0) {
		id = generateCellId();
	}

	if (Cell::getDimension(type) > getDimension()) {
		return cellEnd();
	}

	PiercedVector<Cell>::iterator iterator;
	if (rank == getRank()) {
		iterator = _addInternal(type, std::move(connectStorage), id);
	} else {
		iterator = _addGhost(type, std::move(connectStorage), rank, id);
	}

	return iterator;
}

/*!
	Internal function to add a ghost cell.

	\param type is the type of the cell
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
	\param rank is the rank that owns the cell that will be added
	\param id is the id that will be assigned to the newly created cell.
	If a negative id value is specified, a new unique id will be generated
	for the cell
	\return An iterator pointing to the newly created cell.
*/
PatchKernel::CellIterator PatchKernel::_addGhost(ElementType type, std::unique_ptr<long[]> &&connectStorage,
												 int rank, long id)
{
	// Create the cell
	//
	// If there are internal cells, the ghost cell should be inserted
	// after the last internal cell.
	PiercedVector<Cell>::iterator iterator;
	if (m_lastInternalId < 0) {
		iterator = m_cells.emreclaim(id, id, type, std::move(connectStorage), false, true);
	} else {
		iterator = m_cells.emreclaimAfter(m_lastInternalId, id, id, type, std::move(connectStorage), false, true);
	}
	m_nGhosts++;

	// Update the id of the first ghost cell
	if (m_firstGhostId < 0) {
		m_firstGhostId = id;
	} else if (m_cells.rawIndex(m_firstGhostId) > m_cells.rawIndex(id)) {
		m_firstGhostId = id;
	}

	// Set owner
	setGhostOwner(id, rank);

	return iterator;
}

/*!
	Resore the cell with the specified id.

	The kernel should already contain the cell, only the contents of the
	cell will be updated.

	\param type is the type of the cell
	\param connectivity is the connectivity of the cell
	\param rank is the rank that owns the cell that will be restored
	\param id is the id of the cell that will be restored
	\return An iterator pointing to the restored cell.
*/
PatchKernel::CellIterator PatchKernel::restoreCell(ElementType type, std::unique_ptr<long[]> &&connectStorage,
												   int rank, const long &id)
{
	if (Cell::getDimension(type) > getDimension()) {
		return cellEnd();
	}

	PiercedVector<Cell>::iterator iterator = m_cells.find(id);
	if (iterator == m_cells.end()) {
		throw std::runtime_error("Unable to restore the specified cell: the kernel doesn't contain an entry for that cell.");
	}

	if (rank == getRank()) {
		_restoreInternal(iterator, type, std::move(connectStorage));
	} else {
		_restoreGhost(iterator, type, std::move(connectStorage), rank);
	}

	return iterator;
}

/*!
	Internal function to restore a ghost cell.

	The kernel should already contain the cell, only the contents of the
	cell will be updated.

	\param iterator is an iterator pointing to the cell to restore
	\param type is the type of the cell
	\param connectStorage is the storage the contains or will contain
	the connectivity of the element
	\param rank is the rank that owns the cell that will be restored
*/
void PatchKernel::_restoreGhost(CellIterator iterator, ElementType type,
								std::unique_ptr<long[]> &&connectStorage, int rank)
{
	// Restore cell
	Cell &cell = *iterator;
	cell.initialize(iterator.getId(), type, std::move(connectStorage), false, true);
	m_nGhosts++;

	// Set owner
	setGhostOwner(cell.getId(), rank);
}

/*!
	Internal function to delete a ghost cell.

	\param id is the id of the cell
	\param delayed is true a delayed delete will be performed
*/
void PatchKernel::_deleteGhost(long id, bool delayed)
{
	// Unset ghost owner
	unsetGhostOwner(id);

	// Delete cell
	m_cells.erase(id, delayed);
	m_nGhosts--;
	if (id == m_firstGhostId) {
		updateFirstGhostId();
	}
}

/*!
    Returns iterator to the first ghost cells within the cell list.

    \result An iterator to the first ghost cell.
*/
PatchKernel::CellIterator PatchKernel::ghostBegin()
{
	if (m_nGhosts > 0) {
		return m_cells.find(m_firstGhostId);
	} else {
		return m_cells.end();
	}
}

/*!
	Returns iterator to the end of the list of ghost cells.

	\result An iterator to the end of the list of ghost cell.
*/
PatchKernel::CellIterator PatchKernel::ghostEnd()
{
	return m_cells.end();
}

/*!
    Returns a constant iterator to the first ghost cells within the cell list.

    \result A constant iterator to the first ghost cell.
*/
PatchKernel::CellConstIterator PatchKernel::ghostConstBegin() const
{
	if (m_nGhosts > 0) {
		return m_cells.find(m_firstGhostId);
	} else {
		return m_cells.cend();
	}
}

/*!
	Returns a constant iterator to the end of the list of ghost cells.

	\result A constant iterator to the end of the list of ghost cell.
*/
PatchKernel::CellConstIterator PatchKernel::ghostConstEnd() const
{
	return m_cells.cend();
}

/*!
	Updates the id of the first ghost cell.
*/
void PatchKernel::updateFirstGhostId()
{
	if (m_nGhosts == 0) {
		m_firstGhostId = Cell::NULL_ID;
	} else if (m_nInternals == 0) {
		m_firstGhostId = m_cells.getSizeMarker(m_nInternals, Cell::NULL_ID);
	} else {
		CellIterator first_ghost_iterator = ++m_cells.find(m_lastInternalId);
		m_firstGhostId = first_ghost_iterator->getId();
	}
}

/*!
	Partitions the patch among the processors. Each cell will be assigned
	to a specific processor according to the specified input.

	\param communicator is the communicator that will be used
	\param cellRanks are the ranks of the cells after the partitioning
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\param squeezeStorage if set to true the vector that store patch information
	will be squeezed after the synchronization
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partition(MPI_Comm communicator, const std::vector<int> &cellRanks, bool trackPartitioning, bool squeezeStorage, std::size_t haloSize)
{
	setCommunicator(communicator);

	setHaloSize(haloSize);

	return partition(cellRanks, trackPartitioning, squeezeStorage);
}

/*!
	Partitions the patch among the processors. Each cell will be assigned
	to a specific processor according to the specified input.

	\param cellRanks are the ranks of the cells after the partitioning.
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\param squeezeStorage if set to true the vector that store patch information
	will be squeezed after the synchronization
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partition(const std::vector<int> &cellRanks, bool trackPartitioning, bool squeezeStorage)
{
	std::vector<adaption::Info> partitioningData;

	// Communicator has to be set
	if (!isCommunicatorSet()) {
		throw std::runtime_error ("There is no communicator set for the patch.");
	}

	// Check partitioning status
	PartitioningStatus partitioningStatus = getPartitioningStatus(true);
	if (partitioningStatus == PARTITIONING_UNSUPPORTED) {
		return partitioningData;
	} else if (partitioningStatus != PARTITIONING_CLEAN) {
		throw std::runtime_error ("A partitioning is already in progress.");
	}

	partitioningPrepare(cellRanks, false);

	partitioningData = partitioningAlter(trackPartitioning, squeezeStorage);

	partitioningCleanup();

	return partitioningData;
}

/*!
	Partitions the patch among the processors. The partitioning is done using
	a criteria that tries to balance the load among the processors.

	\param communicator is the communicator that will be used
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\param squeezeStorage if set to true the vector that store patch information
	will be squeezed after the synchronization
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partition(MPI_Comm communicator, bool trackPartitioning, bool squeezeStorage, std::size_t haloSize)
{
	setCommunicator(communicator);

	setHaloSize(haloSize);

	return partition(trackPartitioning, squeezeStorage);
}

/*!
	Partitions the patch among the processors. The partitioning is done using
	a criteria that tries to balance the load among the processors.

	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\param squeezeStorage if set to true the vector that store patch information
	will be squeezed after the synchronization
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partition(bool trackPartitioning, bool squeezeStorage)
{
	std::vector<adaption::Info> partitioningData;

	// Communicator has to be set
	if (!isCommunicatorSet()) {
		throw std::runtime_error ("There is no communicator set for the patch.");
	}

	// Check partitioning status
	PartitioningStatus partitioningStatus = getPartitioningStatus(true);
	if (partitioningStatus == PARTITIONING_UNSUPPORTED) {
		return partitioningData;
	} else if (partitioningStatus != PARTITIONING_CLEAN) {
		throw std::runtime_error ("A partitioning is already in progress.");
	}

	partitioningPrepare(false);

	partitioningData = partitioningAlter(trackPartitioning, squeezeStorage);

	partitioningCleanup();

	return partitioningData;
}

/*!
	Partitions the patch among the processors. Each cell will be assigned
	to a specific processor according to the specified input.

	\param communicator is the communicator that will be used
	\param cellRanks are the ranks of the cells after the partitioning
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partitioningPrepare(MPI_Comm communicator, const std::vector<int> &cellRanks, bool trackPartitioning, std::size_t haloSize)
{
	setCommunicator(communicator);

	setHaloSize(haloSize);

	return partitioningPrepare(cellRanks, trackPartitioning);
}

/*!
	Partitions the patch among the processors. Each cell will be assigned
	to a specific processor according to the specified input.

	\param cellRanks are the ranks of the cells after the partitioning.
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partitioningPrepare(const std::vector<int> &cellRanks, bool trackPartitioning)
{
	std::vector<adaption::Info> partitioningData;

	// Communicator has to be set
	if (!isCommunicatorSet()) {
		throw std::runtime_error ("There is no communicator set for the patch.");
	}

	// Check partitioning status
	PartitioningStatus partitioningStatus = getPartitioningStatus(true);
	if (partitioningStatus == PARTITIONING_UNSUPPORTED) {
		return partitioningData;
	} else if (partitioningStatus != PARTITIONING_CLEAN) {
		throw std::runtime_error ("A partitioning is already in progress.");
	}

	// Build the list of ids to be sent
	auto cellItr = cellBegin();
	for (int k = 0; k < getInternalCount(); ++k) {
		const int rank = cellRanks[k];
		if (rank == getRank()) {
			cellItr++;
			continue;
		}

		m_partitioningLocalSendList[rank].push_back(cellItr->getId());

		cellItr++;
	}

	// Local senders and receivers
	int nLocalExchanges = m_partitioningLocalSendList.size();

	std::vector<int> localSenders;
	std::vector<int> localReceivers;
	localSenders.reserve(nLocalExchanges);
	localReceivers.reserve(nLocalExchanges);
	for (const auto &entry : m_partitioningLocalSendList) {
		localSenders.push_back(getRank());
		localReceivers.push_back(entry.first);
	}

	// Prepare the communication for exchanging the sender/receiver pairs
	std::vector<int> globalExchangeSizes(getProcessorCount());
	MPI_Allgather(&nLocalExchanges, 1, MPI_INT, globalExchangeSizes.data(), 1, MPI_INT, getCommunicator());

	std::vector<int> globalExchangeOffsets(getProcessorCount());
	globalExchangeOffsets[0] = 0;
	for (int i = 1; i < getProcessorCount(); ++i) {
		globalExchangeOffsets[i] = globalExchangeOffsets[i-1] + globalExchangeSizes[i-1];
	}

	// Gather global information
	m_nPartitioningGlobalExchanges = globalExchangeOffsets.back() + globalExchangeSizes.back();

	m_partitioningGlobalSenders.resize(m_nPartitioningGlobalExchanges);
	m_partitioningGlobalReceivers.resize(m_nPartitioningGlobalExchanges);

	MPI_Allgatherv(localSenders.data(), localSenders.size(), MPI_INT, m_partitioningGlobalSenders.data(),
				   globalExchangeSizes.data(), globalExchangeOffsets.data(), MPI_INT, getCommunicator());

	MPI_Allgatherv(localReceivers.data(), localReceivers.size(), MPI_INT, m_partitioningGlobalReceivers.data(),
				   globalExchangeSizes.data(), globalExchangeOffsets.data(), MPI_INT, getCommunicator());

	// Build the information on the cells that will be sent
	if (trackPartitioning) {
		for (const auto &entry : m_partitioningLocalSendList) {
			int receiver = entry.first;
			const std::vector<long> &ids = entry.second;

			partitioningData.emplace_back();
			adaption::Info &partitioningInfo = partitioningData.back();
			partitioningInfo.entity   = adaption::ENTITY_CELL;
			partitioningInfo.type     = adaption::TYPE_PARTITION_SEND;
			partitioningInfo.rank     = receiver;
			partitioningInfo.previous = ids;
		}
	}

	// Update the status
	setPartitioningStatus(PARTITIONING_PREPARED);

	return partitioningData;
}

/*!
	Partitions the patch among the processors. The partitioning is done using
	a criteria that tries to balance the load among the processors.

	\param communicator is the communicator that will be used
	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\param haloSize is the size, expressed in number of layers, of the ghost
	cells halo
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partitioningPrepare(MPI_Comm communicator, bool trackPartitioning, std::size_t haloSize)
{
	setCommunicator(communicator);

	setHaloSize(haloSize);

	return partitioningPrepare(trackPartitioning);
}

/*!
	Partitions the patch among the processors. The partitioning is done using
	a criteria that tries to balance the load among the processors.

	\param trackPartitioning if set to true, the changes to the patch will be
	tracked
	\result Returns a vector of adaption::Info that can be used to track
	the changes done during the partitioning.
*/
std::vector<adaption::Info> PatchKernel::partitioningPrepare(bool trackPartitioning)
{
	std::vector<adaption::Info> partitioningData;

	// Check partitioning status
	PartitioningStatus partitioningStatus = getPartitioningStatus(true);
	if (partitioningStatus == PARTITIONING_UNSUPPORTED) {
		return partitioningData;
	} else if (partitioningStatus != PARTITIONING_CLEAN) {
		throw std::runtime_error ("A partitioning is already in progress.");
	}

	// Reset partitioning information
	m_nPartitioningGlobalExchanges = 0;

	// Execute the partitioning preparation
	partitioningData = _partitioningPrepare(trackPartitioning);

	// Update the status
	setPartitioningStatus(PARTITIONING_PREPARED);

	return partitioningData;
}

/*!
	Alter the patch performing the partitioning.

	The actual modification of the patch takes place during this phase. After
	this phase the adapton is completed and the patch is in its final state.
	Optionally the patch can track the changes performed to the patch.

	\param trackPartitioning if set to true the function will return the changes
	done to the patch during the partitioning
	\param squeezeStorage if set to true patch data structures will be
	squeezed after the partitioning
	\result If the partitioning is tracked, returns a vector of adaption::Info
	with all the changes done to the patch during the partitioning, otherwise
	an empty vector will be returned.
*/
std::vector<adaption::Info> PatchKernel::partitioningAlter(bool trackPartitioning, bool squeezeStorage)
{
	std::vector<adaption::Info> partitioningData;

	// Check partitioning status
	PartitioningStatus partitioningStatus = getPartitioningStatus();
	if (partitioningStatus == PARTITIONING_UNSUPPORTED || partitioningStatus == PARTITIONING_CLEAN) {
		return partitioningData;
	} else if (partitioningStatus != PARTITIONING_PREPARED) {
		throw std::runtime_error ("The prepare function has no been called.");
	}

	// Begin patch alteration
	beginAlteration();

	// Alter patch
	if (m_nPartitioningGlobalExchanges == 0) {
		partitioningData = _partitioningAlter(trackPartitioning);
	} else {
		std::vector<long> emptyCellList;
		for (int i = 0; i < m_nPartitioningGlobalExchanges; ++i) {
			int sender   = m_partitioningGlobalSenders[i];
			int receiver = m_partitioningGlobalReceivers[i];

			std::vector<long> *ids;
			if (sender == getRank()) {
				ids = &(m_partitioningLocalSendList[receiver]);
			} else {
				ids = &emptyCellList;
			}

			adaption::Info partitioningInfo = sendCells_any(sender, receiver, *ids);
			if (trackPartitioning && partitioningInfo.type != adaption::TYPE_NONE) {
				partitioningData.push_back(std::move(partitioningInfo));
			}
		}
	}

	// End patch alteration
	endAlteration(squeezeStorage);

	// The patch is now partitioned
	setPartitioned(true);

	// Update the status
	setPartitioningStatus(PARTITIONING_ALTERED);

	return partitioningData;
}

/*!
	Cleanup patch data structured after the partitioning.

	The patch will only clean-up the data structures needed during the
	partitioning.
*/
void PatchKernel::partitioningCleanup()
{
	PartitioningStatus partitioningStatus = getPartitioningStatus();
	if (partitioningStatus == PARTITIONING_UNSUPPORTED || partitioningStatus == PARTITIONING_CLEAN) {
		return;
	} else if (partitioningStatus == PARTITIONING_PREPARED) {
		throw std::runtime_error ("It is not yet possible to abort a partitioning.");
	} else if (partitioningStatus != PARTITIONING_ALTERED) {
		throw std::runtime_error ("The alter function has no been called.");
	}

	// Clean-up the partitioning
	_partitioningCleanup();

	if (m_nPartitioningGlobalExchanges != 0) {
		m_nPartitioningGlobalExchanges = 0;
		std::unordered_map<int, std::vector<long>>().swap(m_partitioningLocalSendList);
		std::vector<int>().swap(m_partitioningGlobalSenders);
		std::vector<int>().swap(m_partitioningGlobalReceivers);
	}

	// Update the status
	setPartitioningStatus(PARTITIONING_CLEAN);
}

/*!
	Checks if the patch has been partitioned.

	\result Returns true if the patch has been partitioned, false otherwise.
*/
bool PatchKernel::isPartitioned() const
{
	return m_partitioned;
}

/*!
	Sets the partitioned flag.

	\param partitioned is the flag that will be set
*/
void PatchKernel::setPartitioned(bool partitioned)
{
	m_partitioned = partitioned;
}

/*!
	Returns the current partitioning status.

	\param global if set to true the partitioning status will be
	\return The current partitioning status.
*/
PatchKernel::PartitioningStatus PatchKernel::getPartitioningStatus(bool global) const
{
	int partitioningStatus = static_cast<int>(m_partitioningStatus);

	if (global && isCommunicatorSet()) {
		const auto &communicator = getCommunicator();
		MPI_Allreduce(MPI_IN_PLACE, &partitioningStatus, 1, MPI_INT, MPI_MAX, communicator);
	}

	return static_cast<PartitioningStatus>(partitioningStatus);
}

/*!
	Set the current partitioning status.

	\param status is the partitioning status that will be set
*/
void PatchKernel::setPartitioningStatus(PartitioningStatus status)
{
	m_partitioningStatus = status;
}


/*!
	Evaluate partitioning load unbalance index.

	\result Partitioning load unbalance index.
*/
double PatchKernel::evalPartitioningUnbalance()
{
	if (!isPartitioned()) {
		return 0.;
	}

	// Evaluate partition weight
	double localWeight = getInternalCount();

	// Evalaute global weights
	double totalWeight;
	MPI_Allreduce(&localWeight, &totalWeight, 1, MPI_DOUBLE, MPI_SUM, getCommunicator());

	double minimumWeight;
	MPI_Allreduce(&localWeight, &minimumWeight, 1, MPI_DOUBLE, MPI_MIN, getCommunicator());

	double maximumWeight;
	MPI_Allreduce(&localWeight, &maximumWeight, 1, MPI_DOUBLE, MPI_MAX, getCommunicator());

	// Evaluate the unbalance
	double unbalance = (maximumWeight - minimumWeight) / totalWeight;

	return unbalance;
}

/*!
	Prepares the patch for performing the partitioning.

	Default implementation is a no-op function.

	\param trackPartitioning if set to true the function will return the
	changes that will be performed in the alter step
	\result If the partitioning is tracked, returns a vector of adaption::Info
	that can be used to discover what changes will be performed in the alter
	step, otherwise an empty vector will be returned.
*/
std::vector<adaption::Info> PatchKernel::_partitioningPrepare(bool trackPartitioning)
{
	BITPIT_UNUSED(trackPartitioning);

	return std::vector<adaption::Info>();
}

/*!
	Alter the patch performing the partitioning.

	Default implementation is a no-op function.

	\param trackPartitioning if set to true the function will return the changes
	done to the patch during the partitioning
	\result If the partitioning is tracked, returns a vector of adaption::Info
	with all the changes done to the patch during the adaption, otherwise an
	empty vector will be returned.
*/
std::vector<adaption::Info> PatchKernel::_partitioningAlter(bool trackPartitioning)
{
	BITPIT_UNUSED(trackPartitioning);

	assert(false && "The patch needs to implement _partitioningAlter");

	return std::vector<adaption::Info>();
}

/*!
	Cleanup patch data structured after the partitioning.

	Default implementation is a no-op function.
*/
void PatchKernel::_partitioningCleanup()
{
}

/*!
	Gets the rank of the processor that owns the specified cell.

	\param id is the id of the requested cell
	\result The rank that owns the specified cell.
*/
int PatchKernel::getCellRank(const long &id) const
{
	const Cell &cell = getCell(id);
	if (cell.isInterior()) {
		return m_rank;
	} else {
		return m_ghostOwners.at(id);
	}
}

/*!
	Gets the halo layer of the specified cell.

	\param id is the id of the requested cell
	\result The halo layer of the specified cell.
*/
int PatchKernel::getCellHaloLayer(const long &id) const
{
	const Cell &cell = getCell(id);
	if (cell.isInterior()) {
		return 0;
	} else {
		return -1;
	}
}

/*!
	Check if the processors associated to the specified rank is a neighbour.

	\param rank is the rank associated to the processor
	\result True is the processor is a neighbour, false otherwise.
*/
bool PatchKernel::isRankNeighbour(int rank)
{
	return (m_ghostExchangeTargets.count(rank) > 0);
}

/*!
	Get a list of neighbour ranks.

	\result A list of neighbour ranks.
*/
std::vector<int> PatchKernel::getNeighbourRanks()
{
	std::vector<int> neighRanks;
	neighRanks.reserve(m_ghostExchangeTargets.size());
	for (const auto &entry : m_ghostExchangeTargets) {
		neighRanks.push_back(entry.first);
	}

	return neighRanks;
}

/*!
	Gets a constant reference to the ghost targets needed for data exchange.

	\result A constant reference to the ghost targets needed for data
	exchange.
*/
const std::unordered_map<int, std::vector<long>> & PatchKernel::getGhostExchangeTargets() const
{
	return m_ghostExchangeTargets;
}

/*!
	Gets a constant reference to the ghost targets needed for data
	exchange for the specified rank.

	\param rank is the rank for which the information will be retreived
	\result A constant reference to the ghost targets needed for data
	exchange for the specified rank.
*/
const std::vector<long> & PatchKernel::getGhostExchangeTargets(int rank) const
{
	return m_ghostExchangeTargets.at(rank);
}

/*!
	Gets a constant reference to the ghost sources needed for data exchange.

	\result A constant reference to the ghost sources needed for data
	exchange.
*/
const std::unordered_map<int, std::vector<long>> & PatchKernel::getGhostExchangeSources() const
{
	return m_ghostExchangeSources;
}

/*!
	Gets a constant reference to the ghost sources needed for data
	exchange for the specified rank.

	\param rank is the rank for which the information will be retreived
	\result A constant reference to the ghost sources needed for data
	exchange for the specified rank.
*/
const std::vector<long> & PatchKernel::getGhostExchangeSources(int rank) const
{
	return m_ghostExchangeSources.at(rank);
}

/*!
	Sets the owner of the specified ghost.

	\param id is the id of the ghost cell
	\param rank is the rank of the processors that owns the ghost cell
*/
void PatchKernel::setGhostOwner(int id, int rank)
{
	auto ghostOwnerItr = m_ghostOwners.find(id);
	if (ghostOwnerItr != m_ghostOwners.end()) {
		ghostOwnerItr->second = rank;
	} else {
		m_ghostOwners.insert({id, rank});
	}
}

/*!
	Unsets the owner of the specified ghost.

	\param id is the id of the ghost cell
*/
void PatchKernel::unsetGhostOwner(int id)
{
	auto ghostOwnerItr = m_ghostOwners.find(id);
	if (ghostOwnerItr == m_ghostOwners.end()) {
		return;
	}

	m_ghostOwners.erase(ghostOwnerItr);
}

/*!
	Clear the owners of all the ghosts.

	\param updateExchangeInfo if set to true exchange info will be updated
*/
void PatchKernel::clearGhostOwners()
{
	m_ghostOwners.clear();
}

/*!
	Update the information needed for ghost data exchange.
*/
void PatchKernel::updateGhostExchangeInfo()
{
	// Check if all structures needed are ready
	assert(getAdjacenciesBuildStrategy() != ADJACENCIES_NONE);

	// Clear targets
	m_ghostExchangeTargets.clear();

	// Update targets
	for (const auto &entry : m_ghostOwners) {
		int ghostRank = entry.second;
		long ghostId = entry.first;
		m_ghostExchangeTargets[ghostRank].push_back(ghostId);
	}

	// Sort the targets
	for (auto &entry : m_ghostExchangeTargets) {
		std::vector<long> &rankTargets = entry.second;
		std::sort(rankTargets.begin(), rankTargets.end(), CellPositionLess(*this));
	}

	// Clear the sources
	m_ghostExchangeSources.clear();

	// Build the sources
	for (auto &entry : m_ghostExchangeTargets) {
		int rank = entry.first;

		// Generate the source list
		std::vector<long> rankSources = _findGhostExchangeSources(rank);
		if (rankSources.empty()) {
			m_ghostExchangeSources.erase(rank);
			continue;
		}

		// Sort the sources
		std::sort(rankSources.begin(), rankSources.end(), CellPositionLess(*this));

		// Store list
		m_ghostExchangeSources[rank] = std::move(rankSources);
	}
}

/*!
	Finds the internal cells that will be ghost cells for the processors
	with the specified ranks. During data exchange, these cells will be
	the sources form which data will be read from.

	\param rank is the rank for which the information will be built
*/
std::vector<long> PatchKernel::_findGhostExchangeSources(int rank)
{
	// Get targets for the specified rank
	//
	// If there are no targets, there will be no soruces either.
	auto ghostExchangeTargetsItr = m_ghostExchangeTargets.find(rank);
	if (ghostExchangeTargetsItr == m_ghostExchangeTargets.end()) {
		return std::vector<long>();
	}

	std::vector<long> &rankTargets = ghostExchangeTargetsItr->second;

	// The internal neighbours of the ghosts will be sources for the rank
	std::vector<long> neighIds;
	std::unordered_set<long> exchangeSources;
	exchangeSources.reserve(rankTargets.size());
	for (long ghostId : rankTargets) {
		neighIds.clear();
		findCellNeighs(ghostId, &neighIds);
		for (long neighId : neighIds) {
			if (m_ghostOwners.count(neighId) > 0) {
				continue;
			}

			exchangeSources.insert(neighId);
		}
	}

	return std::vector<long>(exchangeSources.begin(), exchangeSources.end());
}

/*!
    Sends the specified list of cells from process with rank sendRank (sender)
    to process with rank recvRank (receiver). If the rank the process currently
    hosting the mesh is neither the sender or the receiver, a notification is
    received in case ghost cells has changed owner.

    \param[in] sendRank sender rank
    \param[in] recvRank receiver rank
    \param[in] cellsToSend list of cells to be moved
    \param[in] squeezeStorage if set to true the vector that store patch information
    will be squeezed after the synchronization
 */
adaption::Info PatchKernel::sendCells(int sendRank, int recvRank,
                                      const std::vector<long> &cellsToSend,
                                      bool squeezeStorage)
{
	//
	// Pereare partitioning alteration
	//
	partitioningPrepare(false);

	//
	// Alter partitioning
	//

	// Begin patch alteration
	beginAlteration();

	// Send cells
	adaption::Info adaptionInfo = sendCells_any(sendRank, recvRank, cellsToSend);

	// End patch alteration
	endAlteration(squeezeStorage);

	// The patch is now partitioned
	setPartitioned(true);

	// Update the status
	setPartitioningStatus(PARTITIONING_ALTERED);

	//
	// Cleanup partitioning alteration
	//
	partitioningCleanup();

	return adaptionInfo;
}

/*!
    Sends the specified list of cells from process with rank sendRank (sender)
    to process with rank recvRank (receiver). If the rank the process currently
    hosting the mesh is neither the sender or the receiver, a notification is
    received in case ghost cells has changed owner.

    \param[in] sendRank sender rank
    \param[in] recvRank receiver rank
    \param[in] cellsToSend list of cells to be moved
 */
adaption::Info PatchKernel::sendCells_any(int sendRank, int recvRank,
                                          const std::vector<long> &cellsToSend)
{
	adaption::Info adaptionInfo;
	if (m_rank == sendRank) {
		adaptionInfo = sendCells_sender(recvRank, cellsToSend);
	} else if (m_rank == recvRank) {
		adaptionInfo = sendCells_receiver(sendRank);
	} else {
		adaptionInfo = sendCells_notified(sendRank, recvRank);
	}

	return adaptionInfo;
}

/*!
    Sends the given list of cells to the process with the specified rank.

    \param[in] recvRank is the receiver rank
    \param[in] cellsToSend is the list of cells to be sent
 */
adaption::Info PatchKernel::sendCells_sender(int recvRank, const std::vector<long> &cellsToSend)
{
    std::vector<long> neighIds;

    //
    // Initialize adaption info
    //
    adaption::Info adaptionInfo;
    adaptionInfo.entity = adaption::ENTITY_CELL;
    adaptionInfo.type   = adaption::TYPE_PARTITION_SEND;
    adaptionInfo.rank   = recvRank;

    //
    // Initialize the list of cells that will be communicated
    //
    // For now the list will contain the cells explicitly marked for sending.
    // Later, a layer of surrounding cells will be added.
    //
    // Only internal cells can be sent.
    std::vector<long> cellsToCommunicate;
    std::unordered_map<long, int> cellRankOnReceiver;

    cellsToCommunicate.reserve(cellsToSend.size());
    cellRankOnReceiver.reserve(cellsToSend.size());
    for (long cellId : cellsToSend) {
		const Cell &cell = m_cells[cellId];
		if (!cell.isInterior()) {
			throw std::runtime_error ("Only internal cells can sent.");
		}

		cellsToCommunicate.push_back(cellId);
        cellRankOnReceiver.insert({{cellId, recvRank}});
    }

    // Identify cells-to-send halo and cells-to-send frame
    //
    // Along with the cells explicitly marked for sending, we need to send an
    // halo of surrounfing cells. Those cells will be used by the receiver to
    // connect the cells it receives to the existing cells, plus some of them
    // will become ghost cells. To identify the cells of the halo, first we
    // find the cells to be send that have at least one face neighbour that
    // remains on this rank and then we find the neighbours of these sends.
    //
    // Cells-to-send frame is made of cells explicitly marked for sending
    // that have at least one neighbour (face, edge, vertex) that will not
    // be sent.

    // Initialize cell frame with cells explicitly marked for sending that
    // have at least one face neighbour that will not be sent
    std::unordered_set<long> cellsToSendFrame;
    for (long cellId : cellsToSend) {
        neighIds.clear();
        findCellFaceNeighs(cellId, &neighIds);
        for (long neighId : neighIds) {
            if (cellRankOnReceiver.count(neighId) == 0) {
                cellsToSendFrame.insert(cellId);
                break;
            }
        }
    }

    // Build cell halo
    std::unordered_set<long> cellsToSendHalo;
    for (long cellId : cellsToSendFrame) {
        neighIds.clear();
        findCellNeighs(cellId, &neighIds);
        for (long neighId : neighIds) {
            if (cellRankOnReceiver.count(neighId) > 0) {
                continue;
            }

            cellsToSendHalo.insert(neighId);
        }
    }

    // At this point, the frame of the cells explicitly marked for sending
    // contains only cells that have one face negihbour that will not be
    // sent. We need to extend the frame to the cells marked for sending
    // that have at least one neighbour (face, edge, vertex) that will not
    // be sent. To find these cells we search among the neighbours of the
    // halo: neighbours that will not be send will be part of the frame.
    for (long cellId : cellsToSendHalo) {
        neighIds.clear();
        findCellNeighs(cellId, &neighIds);
        for (long neighId : neighIds) {
            if (cellRankOnReceiver.count(neighId) > 0) {
                cellsToSendFrame.insert(neighId);
            }
        }
    }

    // Add halo cells to the cells that will be send
    //
    // Cells owned by receiver are already on the receiver, so there is no
    // need to send them.
    for (long cellId : cellsToSendHalo) {
        int ownerRank;
        if (m_ghostOwners.count(cellId) == 0) {
            ownerRank = m_rank;
        } else {
            ownerRank = m_ghostOwners[cellId];
        }

        if (ownerRank != recvRank) {
            cellsToCommunicate.push_back(cellId);
            cellRankOnReceiver.insert({{cellId, ownerRank}});
        }
    }

    //
    // Create the list of vertices to send
    //
    std::unordered_set<long> vertexToCommunicate;
    for (const long cellId : cellsToCommunicate) {
        const Cell &cell = m_cells[cellId];

        ConstProxyVector<long> cellVertexIds = cell.getVertexIds();
        int nCellVertices = cellVertexIds.size();
        for (int j = 0; j < nCellVertices; ++j) {
            long vertexId = cellVertexIds[j];
            if (vertexToCommunicate.count(vertexId) > 0) {
                continue;
            }

            vertexToCommunicate.insert(vertexId);
        }
    }

    //
    // Send vertex data
    //
    OBinaryStream vertexBuffer;
    long vertexBufferSize = 0;

    // Fill buffer with vertex data
    vertexBufferSize += sizeof(long);
    for (long vertexId : vertexToCommunicate) {
        vertexBufferSize += m_vertices[vertexId].getBinarySize();
    }
    vertexBuffer.setSize(vertexBufferSize);

    vertexBuffer << (long) vertexToCommunicate.size();
    for (long vertexId : vertexToCommunicate) {
        vertexBuffer << m_vertices[vertexId];
    }

    if (vertexBufferSize != (long) vertexBuffer.getSize()) {
		throw std::runtime_error ("Cell buffer size does not match calculated size");
	}

    // Communication
    MPI_Send(&vertexBufferSize, 1, MPI_LONG, recvRank, 10, m_communicator);
    MPI_Send(vertexBuffer.data(), vertexBuffer.getSize(), MPI_CHAR, recvRank, 11, m_communicator);

    //
    // Send cell data
    //
    OBinaryStream cellBuffer;
    long cellBufferSize = 0;

    // Fill the buffer with cell data
    cellBufferSize += sizeof(long);
    for (const long cellId : cellsToCommunicate) {
        cellBufferSize += sizeof(int) + sizeof(int) + m_cells[cellId].getBinarySize();
    }
    cellBuffer.setSize(cellBufferSize);

    cellBuffer << (long) cellsToCommunicate.size();
    for (const long cellId : cellsToCommunicate) {
        const Cell &cell = m_cells[cellId];

        // Cells in the cell frame or in the cell halo may already be on
        // receiver. However we don't have enough information to identify
        // those duplicate cells. The receiver needs to perform a check
        // to avoid inserting duplicate cells.
        bool duplicateCheckNeeded = false;
        if (cellsToSendFrame.count(cellId) > 0) {
            duplicateCheckNeeded = true;
        } else if (cellsToSendHalo.count(cellId) > 0) {
            duplicateCheckNeeded = true;
        }
        cellBuffer << duplicateCheckNeeded;

        // Cell owner on receiver
        int cellFutureOwner = cellRankOnReceiver[cellId];
        cellBuffer << cellFutureOwner;

        // Cell data
        cellBuffer << cell;
    }

    if (cellBufferSize != (long) cellBuffer.getSize()) {
		throw std::runtime_error ("Cell buffer size does not match calculated size");
	}

    // Communication
    MPI_Send(&cellBufferSize, 1, MPI_LONG, recvRank, 20, m_communicator);
    MPI_Send(cellBuffer.data(), cellBuffer.getSize(), MPI_CHAR, recvRank, 21, m_communicator);

    //
    // Send notifications
    //
    // Some of the cells that have been send may be ghosts on other processors.
    // Those processors need to be informed of the ownership change of their
    // ghostd.
    //
    // Each processor will send a list of cells that, before the send operation,
    // were owned by the current processor. We need to check which cells have
    // been sent to the sender, i.e., which cells have changed ownership, and
    // notify the ownership change.

    // Build a kd-tree with the centroids of the send cells
    //
    // We need to consider only the cells explicitly marked for sending.
    KdTree<3, Vertex, long> sendCellsCentroidsTree(cellsToSend.size());
    std::vector<Vertex> sendCellsCentroid;
    sendCellsCentroid.reserve(cellsToSend.size());
    for (long cellId : cellsToSend) {
        sendCellsCentroid.emplace_back(Vertex::NULL_ID, evalCellCentroid(cellId));
        sendCellsCentroidsTree.insert(&(sendCellsCentroid.back()), cellId);
    }

    // Receive the number of ownership checks that each processor will request
    int nRanks = getProcessorCount();

    long nOwnershipChecks = 0;
    std::vector<long> nOwnershipChecksPerRank(nRanks);
    MPI_Gather(&nOwnershipChecks, 1, MPI_LONG, nOwnershipChecksPerRank.data(), 1, MPI_LONG, m_rank, m_communicator);

    // Perform ownership checks
    std::vector<int> notifiedRanks;
    for (int rank = 0; rank < nRanks; ++rank) {
        int nRankOwnershipChecks = nOwnershipChecksPerRank[rank];
        if (nRankOwnershipChecks == 0) {
            continue;
        }

        // Receive the centroids that belong to the cells for which the remote
        // rank is requesting an ownership check.
        long recvBufferSize;
        MPI_Recv(&recvBufferSize, 1, MPI_LONG, rank, 30, m_communicator, MPI_STATUS_IGNORE);

        IBinaryStream recvBuffer(recvBufferSize);
        MPI_Recv(recvBuffer.data(), recvBuffer.getSize(), MPI_CHAR, rank, 31, m_communicator, MPI_STATUS_IGNORE);

        // Notify which cells have changed ownership
        //
        // We have received the centroids of cells that, before the send, were
        // owned by this processor. We need to check which ones of these cells
        // have been sent to the receiver.
        long sendBufferSize = nRankOwnershipChecks * sizeof(bool);
        MPI_Send(&sendBufferSize, 1, MPI_LONG, rank, 30, m_communicator);

        OBinaryStream sendBuffer(sendBufferSize);
        for (int k = 0; k < nRankOwnershipChecks; ++k) {
            // Receive remote cell
            Vertex remoteCellCentroid;
            recvBuffer >> remoteCellCentroid;

            // Check if cell has been sent to receiver
            bool cellSentToReceiver = (sendCellsCentroidsTree.exist(&remoteCellCentroid) >= 0);

            // Fill the buffer
            sendBuffer << cellSentToReceiver;
        }
        MPI_Send(sendBuffer.data(), sendBuffer.getSize(), MPI_CHAR, rank, 31, m_communicator);

        // The rank has been notified
        notifiedRanks.push_back(rank);
    }

    //
    // Update cells that no longer belong to this processor
    //

    // Add the sned ids int the adaption info
    //
    // The ids will be sorted by the position of the cells, this is the
    // same order that will be used on the processor that has received the
    // octants. Since the order is the same, the two processors are able
    // to exchange cell data without any additional extra communication
    // (they already know the list of cells for which data is needed and
    // the order in which these data will be sent).
    for (long cellId : cellsToSend) {
        adaptionInfo.previous.push_back(cellId);
    }

    std::sort(adaptionInfo.previous.begin(), adaptionInfo.previous.end(), CellPositionLess(*this));

    // Delete stale cells/interfaces/vertices
    //
    // If the process is senting all its cells we can just clear the patch.
    if (cellsToSend.size() < (std::size_t) getInternalCount()) {
        // Delete sent cells or mark them as ghosts owned by the receiver.
        for (long cellId : cellsToSend) {
            // Check if a cell has to be delete or is a ghost owned by the
            // receiver. A cell will become a ghost if at least one of his
            // neighbours is an internal cell.
            bool moveToGhosts = false;
            if (cellsToSendFrame.count(cellId) > 0) {
                neighIds.clear();
                findCellNeighs(cellId, &neighIds);
                for (long neighId : neighIds) {
                    if (m_ghostOwners.count(neighId) == 0) {
                        moveToGhosts = true;
                        break;
                    }
                }
            }

            // Delete the cell or mark is as a ghost owned by the receiver.
            if (moveToGhosts) {
                moveInternal2Ghost(cellId, recvRank);
            } else {
                deleteCell(cellId, true, true);
            }
        }

        // Delete stale ghosts
        //
        // Loop over all the ghosts and keep only the cells that have at least
        // one internal neighbour.
        auto itr = m_ghostOwners.cbegin();
        while (itr != m_ghostOwners.cend()) {
            long ghostId = itr->first;

            neighIds.clear();
            findCellNeighs(ghostId, &neighIds);
            bool keep = false;
            for (long neighId : neighIds) {
                if (m_ghostOwners.count(neighId) == 0) {
                    keep = true;
                    break;
                }
            }

            auto nextItr = itr;
            nextItr++;
            if (!keep) {
                deleteCell(ghostId, true, true);
            }
            itr = nextItr;
        }

        // Delete orphan interfaces
        deleteOrphanInterfaces();

        // Delete orphan vertices
        deleteOrphanVertices();
    } else {
        // The processor has sent all its cells, the patch is now empty
        reset();
    }

	// Return adaption info
    return adaptionInfo;
}

/*!
    Recevies a list of cells from the specified processor.

    \param[in] sendRank is the rank of the processors sending the cells
 */
adaption::Info PatchKernel::sendCells_receiver(int sendRank)
{
    //
    // Perfom data exchange
    //

    // Vertex data
    long vertexBufferSize;
    MPI_Recv(&vertexBufferSize, 1, MPI_LONG, sendRank, 10, m_communicator, MPI_STATUS_IGNORE);

    IBinaryStream vertexBuffer(vertexBufferSize);
    MPI_Recv(vertexBuffer.data(), vertexBuffer.getSize(), MPI_CHAR, sendRank, 11, m_communicator, MPI_STATUS_IGNORE);

    // Cell data
    long cellBufferSize;
    MPI_Recv(&cellBufferSize, 1, MPI_LONG, sendRank, 20, m_communicator, MPI_STATUS_IGNORE);

    IBinaryStream cellBuffer(cellBufferSize);
    MPI_Recv(cellBuffer.data(), cellBuffer.getSize(), MPI_CHAR, sendRank, 21, m_communicator, MPI_STATUS_IGNORE);

    // Notification data
    long nOwnershipChecks = 0;
    MPI_Gather(&nOwnershipChecks, 1, MPI_LONG, nullptr, 1, MPI_LONG, sendRank, m_communicator);

    //
    // Initialize adaption info
    //
    adaption::Info adaptionInfo;
    adaptionInfo.entity = adaption::ENTITY_CELL;
    adaptionInfo.type   = adaption::TYPE_PARTITION_RECV;
    adaptionInfo.rank   = sendRank;

    //
    // Add vertices
    //

    // Build a kd-tree with the vertices on the ghosts cells
    //
    // These are the only vertices to check for duplicates when receiving the
    // list vertices from the sender.
    //
    // The kd-tree stores the pointer to the vertices. Ifwe try to store in the
    // kd-tree the pointers to the vertices of the patch, the first resize of
    // the vertex container would invalidate the pointer. Create a copy of the
    // vertices and store the pointer to that copy.
    long nGhostsMaxVertices = 0;
    for (const auto &entry : m_ghostOwners) {
        long ghostId = entry.first;
        const Cell &ghost = m_cells[ghostId];
        nGhostsMaxVertices += ghost.getVertexCount();
    }

    KdTree<3, Vertex, long> ghostVerticesTree(nGhostsMaxVertices);
    std::unordered_map<long, Vertex> ghostVertices;
    ghostVertices.reserve(nGhostsMaxVertices);
    for (const auto &entry : m_ghostOwners) {
        long ghostId = entry.first;
        const Cell &ghost = m_cells[ghostId];
        ConstProxyVector<long> ghostVertexIds = ghost.getVertexIds();
        int nGhostVertices = ghostVertexIds.size();
        for (int k = 0; k < nGhostVertices; ++k) {
            long vertexId = ghostVertexIds[k];
            if (ghostVertices.count(vertexId) == 0) {
                ghostVertices.insert({{vertexId, m_vertices[vertexId]}});
                ghostVerticesTree.insert(&ghostVertices.at(vertexId), vertexId);
            }
        }
    }

    // Receive vertices
    //
    // There are no duplicate in the received vertices, but some of them may
    // be already a local vertex of a ghost cell.
    long nRecvVertices;
    vertexBuffer >> nRecvVertices;

    std::unordered_map<long, long> vertexMap;
    vertexMap.reserve(nRecvVertices);
    for (long i = 0; i < nRecvVertices; ++i) {
        Vertex vertex;
        vertexBuffer >> vertex;
        long vertexId = vertex.getId();

        long localVertexId;
        if (ghostVerticesTree.exist(&vertex, localVertexId) < 0) {
        	//Re-use or generate the Id of the vertex
        	if (m_vertexIdGenerator.isAssigned(vertexId)){
        		localVertexId = generateVertexId();
        	}
        	else{
        		localVertexId = vertexId;
        	}
            addVertex(std::move(vertex), localVertexId);
        }
        vertexMap.insert({{vertexId, localVertexId}});
    }

    std::unordered_map<long, Vertex>().swap(ghostVertices);

    //
    // Add cells
    //

    // Receive cells
    //
    // The sender is sending also an halo surroudning the cells explicitly
    // marked for sending. This halo will be used for connecting the received
    // cells to the existing ones and to build the ghosts.
    //
    // Some of the recieved cells may already be among of the ghosts.
    long nReceivedCells;
    cellBuffer >> nReceivedCells;

    std::unordered_map<long, long> cellMap;

    std::vector<long> addedCells;
    addedCells.reserve(nReceivedCells);

    std::vector<long> updateInterfacesCells;
    if (getInterfacesBuildStrategy() == INTERFACES_AUTOMATIC) {
        updateInterfacesCells.reserve(nReceivedCells);
    }

    std::unordered_map<long, FlatVector2D<long>> linkAdjacencies;

    m_cells.reserve(nReceivedCells);
    for (long i = 0; i < nReceivedCells; ++i) {
        // Cell data
        bool duplicateCheckNeeded;
        cellBuffer >> duplicateCheckNeeded;

        int cellOwner;
        cellBuffer >> cellOwner;

        Cell cell;
        cellBuffer >> cell;

        long cellOriginalId = cell.getId();

        // Set cell interior flag
        bool isInterior = (cellOwner == m_rank);
        cell.setInterior(isInterior);

        // Remap connectivity
        cell.renumberVertices(vertexMap);

        // Check if the cells is a duplicate
        //
        // The received cell may be one of the current ghosts.
        long cellId = Cell::NULL_ID;
        if (duplicateCheckNeeded) {
            for (const auto &ghostEntry : m_ghostOwners) {
                const long ghostId = ghostEntry.first;
                const Cell &ghostCell = m_cells[ghostId];
                if (!cell.hasSameConnect(ghostCell)) {
                    continue;
                }

                cellId = ghostId;
                break;
            }
        }

        // If the cell is not a duplicate add it in the cell data structure,
        // otherwise merge the connectivity of the duplicate cell to the
        // existing cell. This ensure that the received cell will be
        // properly connected to the received cells
        if (cellId < 0) {
            // Re-use or generate the Id of the cell
        	if (m_cellIdGenerator.isAssigned(cellOriginalId)){
        		cellId = generateCellId();
        	}
        	else{
        		cellId = cellOriginalId;
        	}

            // Reset the interfaces of the cell, they will be recreated later
            if (getInterfacesBuildStrategy() == INTERFACES_AUTOMATIC) {
                cell.resetInterfaces();
                updateInterfacesCells.push_back(cellId);
            }

            // Add cell
            addCell(std::move(cell), cellOwner, cellId);
            addedCells.push_back(cellId);

            // Update adaption info
            adaptionInfo.current.push_back(cellId);
        } else {
            // Check if the existing cells needs to become an internal cell
            Cell &localCell = m_cells[cellId];
            if (isInterior && !localCell.isInterior()) {
                moveGhost2Internal(cellId);
            }

            // Save the adjacencies of the received cell, this adjacencies
            // will link together the recevied cell to the existing ones.
            FlatVector2D<long> &cellAdjacencies = linkAdjacencies[cellId];

            int nCellFaces = cell.getFaceCount();
            int nCellAdjacencies = cell.getAdjacencyCount();
            cellAdjacencies.reserve(nCellFaces, nCellAdjacencies);
            for (int face = 0; face < nCellFaces; ++face) {
                int nFaceAdjacencies = cell.getAdjacencyCount(face);
                const long *faceAdjacencies = cell.getAdjacencies(face);
                cellAdjacencies.pushBack(nFaceAdjacencies, faceAdjacencies);
            }

            // Delete the interfaces of the cell, they will be recreated later
            if (getInterfacesBuildStrategy() == INTERFACES_AUTOMATIC) {
                // When deleting an interface, the list of cell interfaces is
                // update as well. Therefore, we need to delete the interfaces
                // from the back of the list using a while loop.
                int nLocalCellInterfaces = localCell.getInterfaceCount();
                const long *localCellInterfaces = localCell.getInterfaces();
                while (nLocalCellInterfaces > 0) {
                    const long interfaceId = localCellInterfaces[nLocalCellInterfaces - 1];
                    deleteInterface(interfaceId, true, true);
                    --nLocalCellInterfaces;
                }

                updateInterfacesCells.push_back(cellId);
            }
        }

        // Add the cell to the cell map
        cellMap.insert({{cellOriginalId, cellId}});
    }

    // Remap adjacencies
    for (auto cellId : addedCells) {
        Cell &cell = m_cells[cellId];

        int nCellFaces = cell.getFaceCount();
        for (int face = 0; face < nCellFaces; ++face) {
            int nFaceAdjacencies = cell.getAdjacencyCount(face);
            const long *faceAdjacencies = cell.getAdjacencies(face);
            for (int k = 0; k < nFaceAdjacencies; ++k) {
                long senderAdjacencyId = faceAdjacencies[k];
                if (cellMap.count(senderAdjacencyId) == 0) {
					cell.deleteAdjacency(face, k);
                    continue;
                } else {
					long localAdjacencyId = cellMap.at(senderAdjacencyId);
					cell.setAdjacency(face, k, localAdjacencyId);
				}
            }
        }
    }

    // Link received cells with the current cells
    for (auto &entry : linkAdjacencies) {
        long cellId = entry.first;
        Cell &cell = m_cells[cellId];

        int nCellFaces = cell.getFaceCount();
        FlatVector2D<long> &cellLinkAdjacencies = entry.second;
        for (int face = 0; face < nCellFaces; ++face) {
            int nFaceLinkAdjacencies = cellLinkAdjacencies.getItemCount(face);
            for (int k = 0; k < nFaceLinkAdjacencies; ++k) {
                // We need to updated the adjacencies only if they are cells
                // that have been send.
                long senderAdjacencyId = cellLinkAdjacencies.getItem(face, k);
                if (cellMap.count(senderAdjacencyId) == 0) {
                    continue;
                }

                // If the send cell is already in the adjacency list there is
                // nothing to update.
                long localAdjacencyId  = cellMap.at(senderAdjacencyId);
                if (cell.findAdjacency(face, localAdjacencyId) >= 0) {
                    continue;
                }

                cell.pushAdjacency(face, localAdjacencyId);
            }
        }
    }

    // Update interfaces
    if (getInterfacesBuildStrategy() == INTERFACES_AUTOMATIC) {
        m_interfaces.flush();

        updateInterfaces(updateInterfacesCells);
    }


    // Sort the ids in the adaption info
    //
    // The ids will be sorted by the position of the cells, this is the
    // same order that will be used on the processor that has received the
    // octants. Since the order is the same, the two processors are able
    // to exchange cell data without any additional extra communication
    // (they already know the list of cells for which data is needed and
    // the order in which these data will be sent).
    std::sort(adaptionInfo.current.begin(), adaptionInfo.current.end(), CellPositionLess(*this));

    // Return adaption info
    return adaptionInfo;
}

/*!
    Notifies the current processor of changes in ghost ownership after a
    cell send operation.

    \param[in] sendRank is the rank of the processor sending the cells
    \param[in] recvRank is the rank of the processor receiving the cells
 */
adaption::Info PatchKernel::sendCells_notified(int sendRank, int recvRank)
{
    adaption::Info adaptionInfo;

    // Check if, among the ghosts, there are cells owned by the sender. Those
    // are the cells that may have changed owner.
    std::size_t nGhosts = getGhostCount();

    std::vector<long> ghostsOnSenderRank;
    std::vector<Vertex> ghostCentroidsOnSenderRank;
    ghostsOnSenderRank.reserve(nGhosts);
    ghostCentroidsOnSenderRank.reserve(nGhosts);
    for (const auto &ghostEntry : m_ghostOwners) {
        int ghostOwner = ghostEntry.second;
        if (ghostOwner == sendRank) {
            long ghostId = ghostEntry.first;
            ghostsOnSenderRank.emplace_back(ghostId);
            ghostCentroidsOnSenderRank.emplace_back(Vertex::NULL_ID, evalCellCentroid(ghostEntry.first));
        }
    }

    long nOwnershipChecks = ghostsOnSenderRank.size();
    MPI_Gather(&nOwnershipChecks, 1, MPI_LONG, nullptr, 1, MPI_LONG, sendRank, m_communicator);

    if (nOwnershipChecks == 0) {
        adaptionInfo.type = adaption::TYPE_NONE;
        return adaptionInfo;
    }

    // Send the centroids of the ghosts owned by sender rank
    long sendBufferSize = 0;
    for (const Vertex &centroid : ghostCentroidsOnSenderRank) {
        sendBufferSize += centroid.getBinarySize();
    }
    MPI_Send(&sendBufferSize, 1, MPI_LONG, sendRank, 30, m_communicator);

    OBinaryStream sendBuffer(sendBufferSize);
    for (const Vertex &centroid : ghostCentroidsOnSenderRank) {
        sendBuffer << centroid;
    }
    MPI_Send(sendBuffer.data(), sendBuffer.getSize(), MPI_CHAR, sendRank, 31, m_communicator);

    // Receive ownership change information
    long recvBufferSize;
    MPI_Recv(&recvBufferSize, 1, MPI_LONG, sendRank, 30, m_communicator, MPI_STATUS_IGNORE);

    IBinaryStream recvBuffer(recvBufferSize);
    MPI_Recv(recvBuffer.data(), recvBuffer.getSize(), MPI_CHAR, sendRank, 31, m_communicator, MPI_STATUS_IGNORE);

    // Initialize adaption info
    adaptionInfo.entity = adaption::ENTITY_CELL;
    adaptionInfo.type   = adaption::TYPE_PARTITION_NOTICE;
    adaptionInfo.rank   = sendRank;

    // Apply ownership changes
    for (long ghostId : ghostsOnSenderRank) {
        bool ghostSentToReceiver;
        recvBuffer >> ghostSentToReceiver;

        if (ghostSentToReceiver) {
            setGhostOwner(ghostId, recvRank);
        }
    }

	// Return adaption info
    return adaptionInfo;
}

}

#endif
