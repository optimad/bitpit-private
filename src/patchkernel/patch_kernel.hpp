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

#ifndef __BITPIT_PATCH_KERNEL_HPP__
#define __BITPIT_PATCH_KERNEL_HPP__

#include <cstddef>
#include <deque>
#include <iostream>
#include <memory>
#if BITPIT_ENABLE_MPI==1
#	include <mpi.h>
#endif
#include <set>
#include <string>
#include <vector>
#include <unordered_map>

#include "bitpit_IO.hpp"
#if BITPIT_ENABLE_MPI==1
#	include "bitpit_communications.hpp"
#endif
#include "bitpit_containers.hpp"

#include "adaption.hpp"
#include "cell.hpp"
#include "interface.hpp"
#include "vertex.hpp"

namespace bitpit {

class PatchKernel : public VTKBaseStreamer {

friend class PatchInfo;
friend class PatchNumberingInfo;
friend class PatchManager;

public:
	typedef PiercedVector<Vertex>::iterator VertexIterator;
	typedef PiercedVector<Cell>::iterator CellIterator;
	typedef PiercedVector<Interface>::iterator InterfaceIterator;

	typedef PiercedVector<Vertex>::const_iterator VertexConstIterator;
	typedef PiercedVector<Cell>::const_iterator CellConstIterator;
	typedef PiercedVector<Interface>::const_iterator InterfaceConstIterator;

	typedef PiercedVector<Vertex>::range VertexRange;
	typedef PiercedVector<Cell>::range CellRange;
	typedef PiercedVector<Interface>::range InterfaceRange;

	typedef PiercedVector<Vertex>::const_range VertexConstRange;
	typedef PiercedVector<Cell>::const_range CellConstRange;
	typedef PiercedVector<Interface>::const_range InterfaceConstRange;

	enum WriteTarget {
		WRITE_TARGET_CELLS_ALL
#if BITPIT_ENABLE_MPI
		, WRITE_TARGET_CELLS_INTERNAL
#endif
	};

	/*!
		Functional for comparing the position of two cells.

		The comparison is made with respect to the cell centroid.
	*/
	struct CellPositionLess
	{
		CellPositionLess(const PatchKernel &patch, bool native = true)
			: m_patch(patch), m_native(native)
		{
		}

		virtual ~CellPositionLess() = default;

		bool operator()(long id_1, long id_2) const
		{
			std::array<double, 3> centroid_1;
			std::array<double, 3> centroid_2;
			if (m_native) {
				centroid_1 = m_patch.evalCellCentroid(id_1);
				centroid_2 = m_patch.evalCellCentroid(id_2);
			} else {
				centroid_1 = m_patch.PatchKernel::evalCellCentroid(id_1);
				centroid_2 = m_patch.PatchKernel::evalCellCentroid(id_2);
			}

			for (int k = 0; k < 3; ++k) {
				if (std::abs(centroid_1[k] - centroid_2[k]) <= m_patch.getTol()) {
					continue;
				}

				return centroid_1[k] < centroid_2[k];
			}

			// If we are here the two cell centroids coincide. It's not
			// possible to define an order for the two cells.
			std::ostringstream stream;
			stream << "It was not possible to define an order for cells " << id_1 << " and " << id_2 << ". ";
			stream << "The two cells have the same centroid.";
			throw std::runtime_error (stream.str());
		}

		const PatchKernel &m_patch;
		bool m_native;
	};

	/*!
		Functional for comparing the position of two cells.

		The comparison is made with respect to the cell centroid.
	*/
	struct CellPositionGreater : private CellPositionLess
	{
		CellPositionGreater(const PatchKernel &patch, bool native = true)
			: CellPositionLess(patch, native)
		{
		}

		bool operator()(long id_1, long id_2) const
		{
			return !CellPositionLess::operator()(id_1, id_2);
		}
	};

	/*!
		Functional for comparing the position of two cells.

		WARNING: the function is faster than the comparison based on cell
		centroid but its result is not indepenedent of the vertex order.

		The comparison is made with respect to the position of cell vertices.
	*/
	struct CellFuzzyPositionLess
	{
		CellFuzzyPositionLess(PatchKernel &patch, bool native = true)
			: m_patch(patch), m_native(native)
		{
		}

		virtual ~CellFuzzyPositionLess() = default;

		bool operator()(long id_1, long id_2) const
		{
			// Select the first vertex of the first cell
			ConstProxyVector<long> cellVertexIds_1 = m_patch.getCell(id_1).getVertexIds();

			std::size_t vertexLocalId_1 = 0;
			long vertexId_1 = cellVertexIds_1[vertexLocalId_1];

			// The vertex of the second cell is choosen as the first vertex on
			// that cell not equal to the selected vertex of the first cell.
			ConstProxyVector<long> cellVertexIds_2 = m_patch.getCell(id_2).getVertexIds();
			std::size_t nCellVertices_2 = cellVertexIds_2.size();

			std::size_t vertexLocalId_2 = 0;
			long vertexId_2 = Vertex::NULL_ID;
			while (vertexLocalId_2 <= nCellVertices_2) {
				vertexId_2 = cellVertexIds_2[vertexLocalId_2];
				if (vertexId_1 != vertexId_2) {
					break;
				}

				++vertexLocalId_2;
			}

			// Compare the two vertices
			if (vertexId_2 != Vertex::NULL_ID) {
				const std::array<double, 3> &vertexCoords_1 = m_patch.getVertex(vertexId_1).getCoords();
				const std::array<double, 3> &vertexCoords_2 = m_patch.getVertex(vertexId_2).getCoords();
				for (int k = 0; k < 3; ++k) {
					if (utils::DoubleFloatingEqual()(vertexCoords_1[k], vertexCoords_2[k], m_patch.getTol())) {
						continue;
					}

					return vertexCoords_1[k] < vertexCoords_2[k];
				}
			}

			// If we are here it was not possible to find a vertex on the
			// second cell for the comparison.
			std::ostringstream stream;
			stream << "Unable to fuzzy order cells " << id_1 << " and " << id_2 << ". ";
			throw std::runtime_error (stream.str());
		}

		PatchKernel &m_patch;
		bool m_native;
	};

	/*!
		Functional for comparing the position of two cells.

		WARNING: the function is faster than the comparison based on cell
		centroid but its result is not indepenedent of the vertex order.

		The comparison is made with respect to the position of cell vertices.
	*/
	struct CellFuzzyPositionGreater : private CellFuzzyPositionLess
	{
		CellFuzzyPositionGreater(PatchKernel &patch, bool native = true)
			: CellFuzzyPositionLess(patch, native)
		{
		}

		bool operator()(long id_1, long id_2) const
		{
			return !CellFuzzyPositionLess::operator()(id_1, id_2);
		}
	};

	/*!
		Adjacencies build strategy
	*/
	enum AdjacenciesBuildStrategy {
		ADJACENCIES_NONE = -1,
		ADJACENCIES_AUTOMATIC
	};

	/*!
		Interfaces build strategy
	*/
	enum InterfacesBuildStrategy {
		INTERFACES_NONE = -1,
		INTERFACES_AUTOMATIC
	};

	/*!
		Spawn status
	*/
	enum SpawnStatus {
		SPAWN_UNNEEDED = -1,
		SPAWN_NEEDED,
		SPAWN_DONE
	};

	/*!
		Adaption status
	*/
	enum AdaptionStatus {
		ADAPTION_UNSUPPORTED = -1,
		ADAPTION_CLEAN,
		ADAPTION_DIRTY,
		ADAPTION_PREPARED,
		ADAPTION_ALTERED
	};

	/*!
		Partitioning status
	*/
	enum PartitioningStatus {
		PARTITIONING_UNSUPPORTED = -1,
		PARTITIONING_CLEAN,
		PARTITIONING_PREPARED,
		PARTITIONING_ALTERED
	};

	virtual ~PatchKernel();

	template<typename patch_t>
	static std::unique_ptr<patch_t> clone(const patch_t *original);

	virtual std::unique_ptr<PatchKernel> clone() const = 0;

	virtual void reset();
	virtual void resetVertices();
	virtual void resetCells();
	virtual void resetInterfaces();

	bool reserveVertices(size_t nVertices);
	bool reserveCells(size_t nCells);
	bool reserveInterfaces(size_t nInterfaces);

	std::vector<adaption::Info> update(bool trackAdaption = true, bool squeezeStorage = false);

	virtual void simulateCellUpdate(const long id, adaption::Marker marker, std::vector<Cell> *virtualCells, PiercedVector<Vertex, long> *virtualVertices) const;

	SpawnStatus getSpawnStatus() const;
	std::vector<adaption::Info> spawn(bool trackSpawn);

	AdaptionStatus getAdaptionStatus(bool global = false) const;
	std::vector<adaption::Info> adaption(bool trackAdaption = true, bool squeezeStorage = false);
	std::vector<adaption::Info> adaptionPrepare(bool trackAdaption = true);
	std::vector<adaption::Info> adaptionAlter(bool trackAdaption = true, bool squeezeStorage = false);
	void adaptionCleanup();

	virtual void settleAdaptionMarkers();

	void markCellForRefinement(long id);
	void markCellForCoarsening(long id);
	void resetCellAdaptionMarker(long id);
    adaption::Marker getCellAdaptionMarker(long id);
	void enableCellBalancing(long id, bool enabled);

	bool isDirty(bool global = false) const;
	bool isExpert() const;

	int getId() const;
	int getDimension() const;
	virtual void setDimension(int dimension);
	bool isThreeDimensional() const;

	virtual long getVertexCount() const;
	PiercedVector<Vertex> &getVertices();
	const PiercedVector<Vertex> &getVertices() const;
	Vertex &getVertex(long id);
	const Vertex & getVertex(long id) const;
	const std::array<double, 3> & getVertexCoords(long id) const;
	VertexIterator addVertex(const Vertex &source, long id = Vertex::NULL_ID);
	VertexIterator addVertex(Vertex &&source, long id = Vertex::NULL_ID);
	VertexIterator addVertex(const std::array<double, 3> &coords, long id = Vertex::NULL_ID);
	long countFreeVertices() const;
	long countOrphanVertices() const;
	std::vector<long> findOrphanVertices();
	bool deleteOrphanVertices();
	std::vector<long> collapseCoincidentVertices();
	bool deleteCoincidentVertices();

	VertexIterator getVertexIterator(long id);
	VertexIterator vertexBegin();
	VertexIterator vertexEnd();

	VertexConstIterator getVertexConstIterator(long id) const;
	VertexConstIterator vertexConstBegin() const;
	VertexConstIterator vertexConstEnd() const;

	bool empty() const;

	virtual long getCellCount() const;
	long getInternalCount() const;
#if BITPIT_ENABLE_MPI==1
	long getGhostCount() const;
#endif
	PiercedVector<Cell> &getCells();
	const PiercedVector<Cell> &getCells() const;
	Cell &getCell(long id);
	const Cell &getCell(long id) const;
	virtual ElementType getCellType(long id) const;
	Cell &getLastInternal();
	const Cell &getLastInternal() const;
#if BITPIT_ENABLE_MPI==1
	Cell &getFirstGhost();
	const Cell &getFirstGhost() const;
#endif
	CellIterator addCell(const Cell &source, long id = Element::NULL_ID);
	CellIterator addCell(Cell &&source, long id = Element::NULL_ID);
	CellIterator addCell(ElementType type, long id = Element::NULL_ID);
	CellIterator addCell(ElementType type, const std::vector<long> &connectivity, long id = Element::NULL_ID);
	CellIterator addCell(ElementType type, std::unique_ptr<long[]> &&connectStorage, long id = Element::NULL_ID);
#if BITPIT_ENABLE_MPI==1
	CellIterator addCell(const Cell &source, int rank, long id = Element::NULL_ID);
	CellIterator addCell(Cell &&source, int rank, long id = Element::NULL_ID);
	CellIterator addCell(ElementType type, int rank, long id = Element::NULL_ID);
	CellIterator addCell(ElementType type, const std::vector<long> &connectivity, int rank, long id = Element::NULL_ID);
	CellIterator addCell(ElementType type, std::unique_ptr<long[]> &&connectStorage, int rank, long id = Element::NULL_ID);
#endif
	BITPIT_DEPRECATED(CellIterator addCell(ElementType type, bool interior, long id = Element::NULL_ID));
	BITPIT_DEPRECATED(CellIterator addCell(ElementType type, bool interior, const std::vector<long> &connectivity, long id = Element::NULL_ID));
	BITPIT_DEPRECATED(CellIterator addCell(ElementType type, bool interior, std::unique_ptr<long[]> &&connectStorage, long id = Element::NULL_ID));
	bool deleteCell(long id, bool updateNeighs = true, bool delayed = false);
	bool deleteCells(const std::vector<long> &ids, bool updateNeighs = true, bool delayed = false);
#if BITPIT_ENABLE_MPI==1
	CellIterator moveGhost2Internal(long id);
	CellIterator moveInternal2Ghost(long id, int ownerRank);
#endif
	virtual double evalCellSize(long id) const = 0;
	long countFreeCells() const;
	long countOrphanCells() const;
	virtual std::array<double, 3> evalCellCentroid(long id) const;
	virtual void evalCellBoundingBox(long id, std::array<double,3> *minPoint, std::array<double,3> *maxPoint) const;
	ConstProxyVector<std::array<double, 3>> getCellVertexCoordinates(long id, std::array<double, 3> *staticStorage = nullptr) const;
	std::vector<long> findCellNeighs(long id) const;
	void findCellNeighs(long id, std::vector<long> *neighs) const;
	std::vector<long> findCellNeighs(long id, int codimension, bool complete = true) const;
	void findCellNeighs(long id, int codimension, bool complete, std::vector<long> *neighs) const;
	std::vector<long> findCellFaceNeighs(long id) const;
	void findCellFaceNeighs(long id, std::vector<long> *neighs) const;
	std::vector<long> findCellFaceNeighs(long id, int face) const;
	void findCellFaceNeighs(long id, int face, std::vector<long> *neighs) const;
	std::vector<long> findCellEdgeNeighs(long id, bool complete = true) const;
	void findCellEdgeNeighs(long id, bool complete, std::vector<long> *neighs) const;
	std::vector<long> findCellEdgeNeighs(long id, int edge) const;
	void findCellEdgeNeighs(long id, int edge, std::vector<long> *neighs) const;
	std::vector<long> findCellVertexNeighs(long id, bool complete = true) const;
	void findCellVertexNeighs(long id, bool complete, std::vector<long> *neighs) const;
	std::vector<long> findCellVertexNeighs(long id, int vertex) const;
	void findCellVertexNeighs(long id, int vertex, std::vector<long> *neighs) const;
	std::vector<long> findCellVertexOneRing(long id, int vertex) const;
	void findCellVertexOneRing(long id, int vertex, std::vector<long> *neighs) const;
	bool findFaceNeighCell(long cellId, long neighId, int *cellFace, int *cellAdjacencyId) const;
	std::set<int> getInternalPIDs();
	std::vector<long> getInternalsByPID(int pid);

	CellIterator getCellIterator(long id);
	CellIterator cellBegin();
	CellIterator cellEnd();
	CellIterator internalBegin();
	CellIterator internalEnd();
#if BITPIT_ENABLE_MPI==1
	CellIterator ghostBegin();
	CellIterator ghostEnd();
#endif

	CellConstIterator getCellConstIterator(long id) const;
	CellConstIterator cellConstBegin() const;
	CellConstIterator cellConstEnd() const;
	CellConstIterator internalConstBegin() const;
	CellConstIterator internalConstEnd() const;
#if BITPIT_ENABLE_MPI==1
	CellConstIterator ghostConstBegin() const;
	CellConstIterator ghostConstEnd() const;
#endif

	virtual long getInterfaceCount() const;
	PiercedVector<Interface> &getInterfaces();
	const PiercedVector<Interface> & getInterfaces() const;
	Interface &getInterface(long id);
	const Interface &getInterface(long id) const;
	virtual ElementType getInterfaceType(long id) const;
	InterfaceIterator addInterface(const Interface &source, long id = Element::NULL_ID);
	InterfaceIterator addInterface(Interface &&source, long id = Element::NULL_ID);
	InterfaceIterator addInterface(ElementType type, long id = Element::NULL_ID);
	InterfaceIterator addInterface(ElementType type, const std::vector<long> &connectivity, long id = Element::NULL_ID);
	InterfaceIterator addInterface(ElementType type, std::unique_ptr<long[]> &&connectStorage, long id = Element::NULL_ID);
	bool deleteInterface(long id, bool updateNeighs = true, bool delayed = false);
	bool deleteInterfaces(const std::vector<long> &ids, bool updateNeighs = true, bool delayed = false);
	long countFreeInterfaces() const;
	long countOrphanInterfaces() const;
	std::vector<long> findOrphanInterfaces() const;
	bool deleteOrphanInterfaces();
	bool isInterfaceOrphan(long id) const;
	virtual std::array<double, 3> evalInterfaceCentroid(long id) const;
	virtual void evalInterfaceBoundingBox(long id, std::array<double,3> *minPoint, std::array<double,3> *maxPoint) const;
	ConstProxyVector<std::array<double, 3>> getInterfaceVertexCoordinates(long id, std::array<double, 3> *staticStorage = nullptr) const;

	InterfaceIterator getInterfaceIterator(long id);
	InterfaceIterator interfaceBegin();
	InterfaceIterator interfaceEnd();

	InterfaceConstIterator getInterfaceConstIterator(long id) const;
	InterfaceConstIterator interfaceConstBegin() const;
	InterfaceConstIterator interfaceConstEnd() const;

	long countFaces() const;
	long countFreeFaces() const;

	bool sort();
	bool sortVertices();
	bool sortCells();
	bool sortInterfaces();

	bool squeeze();
	bool squeezeVertices();
	bool squeezeCells();
	bool squeezeInterfaces();

	long locatePoint(double x, double y, double z);
	virtual long locatePoint(const std::array<double, 3> &point) = 0;
	bool isSameFace(long cellId_A, int face_A, long cellId_B, int face_B);

	AdjacenciesBuildStrategy getAdjacenciesBuildStrategy() const;
	void clearAdjacencies();
	virtual void buildAdjacencies();
	virtual void updateAdjacencies(const std::vector<long> &cellIds);

	InterfacesBuildStrategy getInterfacesBuildStrategy() const;
	void clearInterfaces();
	virtual void buildInterfaces();
	virtual void updateInterfaces(const std::vector<long> &cellIds);

	void getBoundingBox(std::array<double, 3> &minPoint, std::array<double, 3> &maxPoint) const;
	void getBoundingBox(bool global, std::array<double, 3> &minPoint, std::array<double, 3> &maxPoint) const;
	bool isBoundingBoxDirty(bool global = false) const;
	void updateBoundingBox(bool forcedUpdated = false);

	virtual void translate(std::array<double, 3> translation);
	void translate(double sx, double sy, double sz);
	virtual void scale(std::array<double, 3> scaling);
	void scale(double scaling);
	void scale(double sx, double sy, double sz);

	void setTol(double tolerance);
	double getTol() const;
	void resetTol();
	bool isTolCustomized() const;

	void extractEnvelope(PatchKernel &envelope) const;

	void displayTopologyStats(std::ostream &out, unsigned int padding = 0) const;
	void displayVertices(std::ostream &out, unsigned int padding = 0) const;
	void displayCells(std::ostream &out, unsigned int padding = 0) const;
	void displayInterfaces(std::ostream &out, unsigned int padding = 0) const;

	VTKUnstructuredGrid & getVTK();
	WriteTarget getVTKWriteTarget() const;
	void setVTKWriteTarget(WriteTarget targetCells);
	const PatchKernel::CellConstRange getVTKCellWriteRange() const;
	void write(VTKWriteMode mode = VTKWriteMode::DEFAULT);
	void write(const std::string &name, VTKWriteMode mode = VTKWriteMode::DEFAULT);

	void flushData(std::fstream &stream, const std::string &name, VTKFormat format) override;

	int getDumpVersion() const;
	void dump(std::ostream &stream) const;
	void restore(std::istream &stream, bool reregister = false);

	void consecutiveRenumberVertices(long offset = 0);
	void consecutiveRenumberCells(long offset = 0);
	void consecutiveRenumberInterfaces(long offset = 0);
	void consecutiveRenumber(long offsetVertices, long offsetCells, long offsetInterfaces);

#if BITPIT_ENABLE_MPI==1
	virtual void setCommunicator(MPI_Comm communicator);
	void freeCommunicator();
	bool isCommunicatorSet() const;
	const MPI_Comm & getCommunicator() const;
	int getRank() const;
	int getProcessorCount() const;

	void setHaloSize(std::size_t haloSize);
	std::size_t getHaloSize() const;

	int getCellRank(long id) const;
	virtual int getCellHaloLayer(long id) const;

	bool isRankNeighbour(int rank);
	std::vector<int> getNeighbourRanks();
	const std::unordered_map<int, std::vector<long>> & getGhostExchangeTargets() const;
	const std::vector<long> & getGhostExchangeTargets(int rank) const;
	const std::unordered_map<int, std::vector<long>> & getGhostExchangeSources() const;
	const std::vector<long> & getGhostExchangeSources(int rank) const;

	bool isPartitioned() const;
	PartitioningStatus getPartitioningStatus(bool global = false) const;
	double evalPartitioningUnbalance();
	std::vector<adaption::Info> partition(MPI_Comm communicator, const std::vector<int> &cellRanks, bool trackPartitioning, bool squeezeStorage = false, std::size_t haloSize = 1);
	std::vector<adaption::Info> partition(const std::vector<int> &cellRanks, bool trackPartitioning, bool squeezeStorage = false);
	std::vector<adaption::Info> partition(MPI_Comm communicator, bool trackPartitioning, bool squeezeStorage = false, std::size_t haloSize = 1);
	std::vector<adaption::Info> partition(bool trackPartitioning, bool squeezeStorage = false);
	std::vector<adaption::Info> partitioningPrepare(MPI_Comm communicator, const std::vector<int> &cellRanks, bool trackPartitioning, std::size_t haloSize = 1);
	std::vector<adaption::Info> partitioningPrepare(const std::vector<int> &cellRanks, bool trackPartitioning);
	std::vector<adaption::Info> partitioningPrepare(MPI_Comm communicator, bool trackPartitioning, std::size_t haloSize = 1);
	std::vector<adaption::Info> partitioningPrepare(bool trackPartitioning);
	std::vector<adaption::Info> partitioningAlter(bool trackPartitioning = true, bool squeezeStorage = false);
	void partitioningCleanup();

	adaption::Info sendCells(int sendRank, int recvRank, const std::vector<long> &cellsToSend, bool squeezeStorage = false);
#endif

protected:
	PiercedVector<Vertex> m_vertices;
	PiercedVector<Cell> m_cells;
	PiercedVector<Interface> m_interfaces;

	PatchKernel(bool expert);
	PatchKernel(int dimension, bool expert);
	PatchKernel(int id, int dimension, bool expert);
	PatchKernel(const PatchKernel &other);
    PatchKernel & operator=(const PatchKernel &other) = delete;

	void clearBoundingBox();
	bool isBoundingBoxFrozen() const;
	void setBoundingBoxFrozen(bool frozen);
	void setBoundingBoxDirty(bool dirty);
	void setBoundingBox(const std::array<double, 3> &minPoint, const std::array<double, 3> &maxPoint);

	CellIterator restoreCell(ElementType type, std::unique_ptr<long[]> &&connectStorage, long id);
#if BITPIT_ENABLE_MPI==1
	CellIterator restoreCell(ElementType type, std::unique_ptr<long[]> &&connectStorage, int rank, long id);
#endif

	InterfaceIterator restoreInterface(ElementType type, std::unique_ptr<long[]> &&connectStorage, long id);

	VertexIterator restoreVertex(const std::array<double, 3> &&coords, long id);

	bool deleteVertex(long id, bool delayed = false);
	bool deleteVertices(const std::vector<long> &ids, bool delayed = false);

	void dumpVertices(std::ostream &stream) const;
	void restoreVertices(std::istream &stream);

	void dumpCells(std::ostream &stream) const;
	void restoreCells(std::istream &stream);

	void dumpInterfaces(std::ostream &stream) const;
	void restoreInterfaces(std::istream &stream);

	void updateLastInternalId();
#if BITPIT_ENABLE_MPI==1
	void updateFirstGhostId();
#endif

	std::unordered_map<long, std::vector<long>> binGroupVertices(const PiercedVector<Vertex> &vertices, int nBins);
	std::unordered_map<long, std::vector<long>> binGroupVertices(int nBins);

	void setAdjacenciesBuildStrategy(AdjacenciesBuildStrategy status);

	void setInterfacesBuildStrategy(InterfacesBuildStrategy status);

	void setSpawnStatus(SpawnStatus status);
	virtual std::vector<adaption::Info> _spawn(bool trackAdaption);

	void setAdaptionStatus(AdaptionStatus status);
	virtual std::vector<adaption::Info> _adaptionPrepare(bool trackAdaption);
	virtual std::vector<adaption::Info> _adaptionAlter(bool trackAdaption);
	virtual void _adaptionCleanup();
	virtual bool _markCellForRefinement(long id);
	virtual bool _markCellForCoarsening(long id);
	virtual bool _resetCellAdaptionMarker(long id);
	virtual adaption::Marker _getCellAdaptionMarker(long id);
	virtual bool _enableCellBalancing(long id, bool enabled);

	virtual void _setTol(double tolerance);
	virtual void _resetTol();

	virtual int _getDumpVersion() const = 0;
	virtual void _dump(std::ostream &stream) const = 0;
	virtual void _restore(std::istream &stream) = 0;

	virtual long _getCellNativeIndex(long id) const;

	virtual void _findCellNeighs(long id, const std::vector<long> &blackList, std::vector<long> *neighs) const;
	virtual void _findCellFaceNeighs(long id, int face, const std::vector<long> &blackList, std::vector<long> *neighs) const;
	virtual void _findCellEdgeNeighs(long id, int edge, const std::vector<long> &blackList, std::vector<long> *neighs) const;
	virtual void _findCellVertexNeighs(long id, int vertex, const std::vector<long> &blackList, std::vector<long> *neighs) const;

	void setExpert(bool expert);

	void addPointToBoundingBox(const std::array<double, 3> &point);
	void removePointFromBoundingBox(const std::array<double, 3> &point, bool delayedBoxUpdate = false);
#if BITPIT_ENABLE_MPI==1
	virtual std::size_t _getMaxHaloSize();
	virtual void _setHaloSize(std::size_t haloSize);

	void setPartitioned(bool partitioned);
	void setPartitioningStatus(PartitioningStatus status);
	virtual std::vector<adaption::Info> _partitioningPrepare(bool trackPartitioning);
	virtual std::vector<adaption::Info> _partitioningAlter(bool trackPartitioning);
	virtual void _partitioningCleanup();

	virtual std::vector<long> _findGhostExchangeSources(int rank);
#endif

	template<typename item_t, typename id_t = long>
	std::unordered_map<id_t, id_t> consecutiveItemRenumbering(PiercedVector<item_t, id_t> &container, long offset);

	template<typename item_t, typename id_t = long>
	void mappedItemRenumbering(PiercedVector<item_t, id_t> &container, const std::unordered_map<id_t, id_t> &renumberMap);

	ConstProxyVector<std::array<double, 3>> getElementVertexCoordinates(const Element &element, std::array<double, 3> *externalStorage = nullptr) const;

private:
	double DEFAULT_TOLERANCE = 1e-14;

	IndexGenerator<long> m_vertexIdGenerator;
	IndexGenerator<long> m_interfaceIdGenerator;
	IndexGenerator<long> m_cellIdGenerator;

	long m_nInternals;
#if BITPIT_ENABLE_MPI==1
	long m_nGhosts;
#endif

	long m_lastInternalId;
#if BITPIT_ENABLE_MPI==1
	long m_firstGhostId;
#endif

	VTKUnstructuredGrid m_vtk ;
	WriteTarget m_vtkWriteTarget;
	PiercedStorage<long, long> m_vtkVertexMap;

	bool m_boxFrozen;
	bool m_boxDirty;
	std::array<double, 3> m_boxMinPoint;
	std::array<double, 3> m_boxMaxPoint;
	std::array<int, 3> m_boxMinCounter;
	std::array<int, 3> m_boxMaxCounter;

	AdjacenciesBuildStrategy m_adjacenciesBuildStrategy;

	InterfacesBuildStrategy m_interfacesBuildStrategy;

	SpawnStatus m_spawnStatus;

	AdaptionStatus m_adaptionStatus;

	bool m_expert;

	int m_id;
	int m_dimension;

	bool m_hasCustomTolerance;
	double m_tolerance;

	int m_rank;
	int m_nProcessors;
#if BITPIT_ENABLE_MPI==1
	MPI_Comm m_communicator;
	bool m_partitioned;
	PartitioningStatus m_partitioningStatus;

	std::size_t m_haloSize;

	int m_nPartitioningGlobalExchanges;
	std::vector<int> m_partitioningGlobalSenders;
	std::vector<int> m_partitioningGlobalReceivers;
	std::unordered_map<int, std::vector<long>> m_partitioningLocalSendList;

	std::unordered_map<long, int> m_ghostOwners;
	std::unordered_map<int, std::vector<long>> m_ghostExchangeTargets;
	std::unordered_map<int, std::vector<long>> m_ghostExchangeSources;

	void setGhostOwner(int id, int rank);
	void unsetGhostOwner(int id);
	void clearGhostOwners();

	adaption::Info sendCells_any(int sendRank, int recvRank, const std::vector<long> &cellsToSend);
	adaption::Info sendCells_sender(int recvRank, const std::vector<long> &cellsToSend);
	adaption::Info sendCells_receiver(int sendRank);
	adaption::Info sendCells_notified(int sendRank, int recvRank);

	void updateGhostExchangeInfo();
#endif

	void initialize();

	void beginAlteration();
	void endAlteration(bool squeezeStorage = false);

	InterfaceIterator buildCellInterface(Cell *cell_1, int face_1, Cell *cell_2, int face_2, long interfaceId = Element::NULL_ID);

	int findAdjoinNeighFace(long cellId, long neighId) const;

	void setId(int id);

	std::array<double, 3> evalElementCentroid(const Element &element) const;
	void evalElementBoundingBox(const Element &element, std::array<double,3> *minPoint, std::array<double,3> *maxPoint) const;

	void mergeAdaptionInfo(std::vector<adaption::Info> &&source, std::vector<adaption::Info> &destination);

	CellIterator _addInternal(ElementType type, std::unique_ptr<long[]> &&connectStorage, long id);
#if BITPIT_ENABLE_MPI==1
	CellIterator _addGhost(ElementType type, std::unique_ptr<long[]> &&connectStorage, int rank, long id);
#endif

	void _restoreInternal(CellIterator iterator, ElementType type, std::unique_ptr<long[]> &&connectStorage);
#if BITPIT_ENABLE_MPI==1
	void _restoreGhost(CellIterator iterator, ElementType type, std::unique_ptr<long[]> &&connectStorage, int rank);
#endif

	void _deleteInternal(long id, bool delayed);
#if BITPIT_ENABLE_MPI==1
	void _deleteGhost(long id, bool delayed);
#endif
};

}

// Template implementation
#include "patch_kernel.tpp"

#endif
