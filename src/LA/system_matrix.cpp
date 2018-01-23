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

#include "system_matrix.hpp"

namespace bitpit {

/**
* \class SparseMatrix
* \ingroup system_solver
*
* \brief Sparse matrix.
*/

#if BITPIT_ENABLE_MPI==1
/**
* Default constructor
*
* \param communicator is the MPI communicator
*/
SparseMatrix::SparseMatrix(MPI_Comm communicator)
    : m_communicator(communicator),
      m_maxRowNZ(0), m_lastRow(-1),
      m_global_nRows(0), m_global_nCols(0), m_global_nNZ(0),
      m_global_maxRowNZ(0), m_global_rowOffset(0), m_global_colOffset(0)
{
}
#else
/**
* Default constructor
*/
SparseMatrix::SparseMatrix()
    : m_maxRowNZ(0), m_lastRow(-1)
{
}
#endif

/**
* Constructor
*
*/
#if BITPIT_ENABLE_MPI==1
/*!
* \param communicator is the MPI communicator
*/
#endif
/*!
* \param nRows is the number of rows of the matrix
* \param nCols is the number of columns of the matrix
* \param nMaximumNZ is the maximum number of non-zero elements that the matrix
* will contain. This is just an optional hint. If the actual number of non-zero
* elements turns out to be greater than the provided value, the initialization
* of the matrix pattern will be slower because reallocation of internal data
* may be needed
*/
#if BITPIT_ENABLE_MPI==1
SparseMatrix::SparseMatrix(MPI_Comm communicator, long nRows, long nCols, long nMaximumNZ)
    : SparseMatrix(communicator)
#else
SparseMatrix::SparseMatrix(long nRows, long nCols, long nMaximumNZ)
    : SparseMatrix()
#endif
{
    _initialize(nRows, nCols, nMaximumNZ);
}

/**
* Initialize the pattern.
*
* \param nRows is the number of local rows of the matrix
* \param nCols is the number of local columns of the matrix
* \param nMaximumNZ is the maximum number of non-zero elements that the matrix
* will contain. This is just an optional hint. If the actual number of non-zero
* elements turns out to be greater than the provided value, the initialization
* of the matrix pattern will be slower because reallocation of internal data
* may be needed
*/
void SparseMatrix::initialize(long nRows, long nCols, long nMaximumNZ)
{
    _initialize(nRows, nCols, nMaximumNZ);
}

/**
* Internal function to initialize the pattern.
*
* \param nRows is the number of local rows of the matrix
* \param nCols is the number of local columns of the matrix
* \param nMaximumNZ is the maximum number of non-zero elements that the matrix
* will contain. This is just an optional hint. If the actual number of non-zero
* elements turns out to be greater than the provided value, the initialization
* of the matrix pattern will be slower because reallocation of internal data
* may be needed
*/
void SparseMatrix::_initialize(long nRows, long nCols, long nMaximumNZ)
{
    assert(nRows >= 0);
    assert(nCols >= 0);
    assert(nMaximumNZ >= 0);

    clear();

    m_pattern.reserve(nRows, nMaximumNZ);
    m_values.reserve(nMaximumNZ);

    m_nRows = nRows;
    m_nCols = nCols;
#if BITPIT_ENABLE_MPI == 1
    MPI_Allreduce(&m_nRows, &m_global_nRows, 1, MPI_LONG, MPI_SUM, m_communicator);
    MPI_Allreduce(&m_nCols, &m_global_nCols, 1, MPI_LONG, MPI_SUM, m_communicator);

    int nProcessors;
    MPI_Comm_size(m_communicator, &nProcessors);

    m_global_rowOffset = 0;
    m_global_colOffset = 0;
    if (nProcessors > 1) {
        std::vector<long> nGlobalRows(nProcessors);
        MPI_Allgather(&nRows, 1, MPI_LONG, nGlobalRows.data(), 1, MPI_LONG, m_communicator);

        std::vector<long> nGlobalCols(nProcessors);
        MPI_Allgather(&nCols, 1, MPI_LONG, nGlobalCols.data(), 1, MPI_LONG, m_communicator);

        int rank;
        MPI_Comm_rank(m_communicator, &rank);
        for (int i = 0; i < rank; ++i) {
            m_global_rowOffset += nGlobalRows[i];
            m_global_colOffset += nGlobalCols[i];
        }
    }
#else
    BITPIT_UNUSED(nCols)
#endif
}


/**
* Clear the pattern.
*
* \param release if set to true the memory hold by the pattern will be released
*/
void SparseMatrix::clear(bool release)
{
    if (release) {
        m_pattern.clear();
        m_values.clear();
    } else {
        FlatVector2D<long>().swap(m_pattern);
        std::vector<double>().swap(m_values);
    }

    m_nRows    = -1;
    m_nCols    = -1;
    m_nNZ      = -1;
    m_maxRowNZ = -1;
    m_lastRow  = -1;

#if BITPIT_ENABLE_MPI==1
    m_global_nRows     = 0;
    m_global_nCols     = 0;
    m_global_nNZ       = 0;
    m_global_maxRowNZ  = 0;
    m_global_rowOffset = 0;
    m_global_colOffset = 0;
#endif
}

/*!
* Squeeze.
*
* Requests the matrix pattern to reduce its capacity to fit its size.
*/
void SparseMatrix::squeeze()
{
    m_pattern.shrinkToFit();
    m_values.shrink_to_fit();
}

/**
* Check if all the rows of the pattern are defined. If this is the case, the
* pattern is considered frozen and no more rows can be added.

* \result Returns true if the pattern is finalized, false otherwise.
*/
bool SparseMatrix::isInitialized() const
{
    return (m_nRows >= 0);
}

/**
* Check if all the rows of the pattern are defined. If this is the case, the
* pattern is considered frozen and no more rows can be added.

* \result Returns true if the pattern is finalized, false otherwise.
*/
bool SparseMatrix::isFinalized() const
{
    std::cout << " SS " << m_lastRow  << std::endl;
    std::cout << " SS " << m_nRows << std::endl;
    return (m_lastRow >= (m_nRows - 1));
}

/**
* Get the number of rows of the matrix.
*
* \result The number of rows of the matrix.
*/
long SparseMatrix::getRowCount() const
{
    return m_nRows;
}

/**
* Get the number of columns of the matrix.
*
* \result The number of columns of the matrix.
*/
long SparseMatrix::getColCount() const
{
    return m_nCols;
}

/**
* Get the number of non-zero elements.
*
* \result The number of non-zero elements.
*/
long SparseMatrix::getNZCount() const
{
    return m_nNZ;
}

/**
* Get the number of non-zero elements in the specified row.
*
* \result The number of non-zero elements in the specified row.
*/
long SparseMatrix::getRowNZCount(long row) const
{
    return m_pattern.getItemCount(row);
}

/**
* Get the maximum number of non-zero elements per row.
*
* \result The maximum number of non-zero elements per row.
*/
long SparseMatrix::getMaxRowNZCount() const
{
    return m_maxRowNZ;
}

#if BITPIT_ENABLE_MPI==1
/**
* Get the number of global rows
*
* \result The number of global rows
*/
long SparseMatrix::getRowGlobalCount() const
{
    return m_global_nRows;
}

/**
* Get number of global columns.
*
* \result The number of global columns.
*/
long SparseMatrix::getColGlobalCount() const
{
    return m_global_nCols;
}

/**
* Get the global number of non-zero elements.
*
* \result The global number of non-zero elements.
*/
long SparseMatrix::getNZGlobalCount() const
{
    return m_global_nNZ;
}

/**
* Get the global maximum number of non-zero elements per row.
*
* \result The global maximum number of non-zero elements per row.
*/
long SparseMatrix::getMaxRowNZGlobalCount() const
{
    return m_global_maxRowNZ;
}

/**
* Extract the list of local global rows.
*
* \result The the list of local global rows.
*/
std::vector<long> SparseMatrix::extractLocalGlobalRows() const
{
    const std::size_t *rowExtents = m_pattern.indices();

    std::vector<long> localGlobalRows;
    localGlobalRows.reserve(m_nRows);
    for (long i = 0; i < m_nRows; ++i) {
        std::size_t nRowNZ = rowExtents[i + 1] - rowExtents[i];
        if (nRowNZ > 0) {
            localGlobalRows.push_back(m_global_rowOffset + i);
        }
    }

    return localGlobalRows;
}

/**
* Extract the list of ghost global rows.
*
* \result The the list of ghost global rows.
*/
std::vector<long> SparseMatrix::extractGhostGlobalRows() const
{
    return std::vector<long>();
}

/**
* Extract the list of local global columns.
*
* \result The the list of local global columns.
*/
std::vector<long> SparseMatrix::extractLocalGlobalCols() const
{
    long firstGlobalCol = m_global_colOffset;
    long lastGlobalCol  = m_global_colOffset + m_nCols - 1;

    const long *globalCols = m_pattern.data();

    std::vector<long> localGlobalCols;
    localGlobalCols.reserve(m_nNZ);
    for (long k = 0; k < m_nNZ; ++k) {
        long globalCol = globalCols[k];
        if (globalCol < firstGlobalCol || globalCol >= lastGlobalCol) {
            utils::addToOrderedVector<long>(globalCol, localGlobalCols);
        }
    }

    return localGlobalCols;
}

/**
* Extract the list of ghost global columns.
*
* \result The the list of ghost global columns.
*/
std::vector<long> SparseMatrix::extractGhostGlobalCols() const
{
    long firstGlobalCol = m_global_colOffset;
    long lastGlobalCol  = m_global_colOffset + m_nCols - 1;

    const long *globalCols = m_pattern.data();

    std::vector<long> ghostGlobalCols;
    ghostGlobalCols.reserve(m_nNZ);
    for (long k = 0; k < m_nNZ; ++k) {
        long globalCol = globalCols[k];
        if (globalCol >= firstGlobalCol || globalCol < lastGlobalCol) {
            utils::addToOrderedVector<long>(globalCol, ghostGlobalCols);
        }
    }

    return ghostGlobalCols;
}
#endif

/**
* Add a row.
*
* \param rowPattern are the indexes of the non-zero columns of the matrix
* \param rowValues are the values of the non-zero columns of the matrix
*/
void SparseMatrix::addRow(const std::vector<long> &rowPattern, const std::vector<double> &rowValues)
{
    addRow(rowPattern.size(), rowPattern.data(), rowValues.data());
}

/**
* Add a row.
*
* \param nRowNZ is the number of non-zero elements in the row
* \param rowPattern are the indexes of the non-zero columns of the matrix
* \param rowValues are the values of the non-zero columns of the matrix
*/
void SparseMatrix::addRow(long nRowNZ, const long *rowPattern, const double *rowValues)
{
    if (isFinalized()) {
        throw std::runtime_error("Unable to add another row: all rows have already been defined.");
    }

    // Add the row pattern
    m_pattern.pushBack(nRowNZ, rowPattern);

    // Add row values
    m_values.insert(m_values.end(), rowValues, rowValues + nRowNZ);

    // Update the non-zero element conuters
    m_maxRowNZ = std::max(nRowNZ, m_maxRowNZ);
    m_nNZ += nRowNZ;

    // Update the index of the last row
    //
    // When inserting the last row the matrix can be finalized.
    m_lastRow++;
    if (isFinalized()) {
#if BITPIT_ENABLE_MPI==1
        // Updathe global information of the non-zero elements
        MPI_Allreduce(&m_maxRowNZ, &m_global_maxRowNZ, 1, MPI_LONG, MPI_MAX, m_communicator);
        MPI_Allreduce(&m_nNZ, &m_global_nNZ, 1, MPI_LONG, MPI_SUM, m_communicator);
#endif
    }
}

/**
* Get the pattern of the specified row.
*
* \param row is the row
* \result The pattern of the specified row.
*/
ConstProxyVector<long> SparseMatrix::getRowPattern(long row) const
{
    const std::size_t *rowExtent = m_pattern.indices(row);

    const long *rowPattern = m_pattern.data() + rowExtent[0];
    const std::size_t rowPatternSize = rowExtent[1] - rowExtent[0];

    return ConstProxyVector<long>(rowPattern, rowPatternSize);
}

/**
* Get the values of the specified row.
*
* \param row is the row
* \result The values of the specified row.
*/
ConstProxyVector<double> SparseMatrix::getRowValues(long row) const
{
    const std::size_t *rowExtent = m_pattern.indices(row);

    const double *rowValues = m_values.data() + rowExtent[0];
    const std::size_t nRowValues = rowExtent[1] - rowExtent[0];

    return ConstProxyVector<double>(rowValues, nRowValues);
}

}
