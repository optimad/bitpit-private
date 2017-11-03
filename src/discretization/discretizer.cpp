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

#include <stdexcept>
#include <string>

#include "discretizer.hpp"

namespace bitpit {

/*!
 * Default constuctor
 */
#if BITPIT_ENABLE_MPI==1
Discretizer::Discretizer(MPI_Comm communicator, bool debug)
    : SystemSolver(communicator, debug), m_matrix(communicator)
#else
Discretizer::Discretizer(bool debug)
    : SystemSolver(debug), m_matrix
#endif
{
}

/*!
 * Constuctor
 */
#if BITPIT_ENABLE_MPI==1
Discretizer::Discretizer(MPI_Comm communicator, long nUnknowns, long nMaximumNZ, bool debug)
    : SystemSolver(communicator, debug), m_matrix(communicator, nUnknowns, nUnknowns, nMaximumNZ)
#else
Discretizer::Discretizer(long nUnknowns, long nMaximumNZ, bool debug)
    : SystemSolver(debug), m_matrix(nUnknowns, nUnknowns, nMaximumNZ)
#endif
{
    _initialize(nUnknowns, nMaximumNZ);
}

/*!
* Initialize the discretizer.
*
* \param nUnknowns is the number of unknowns
* \param nMaximumNZ is the maximum number of non-zero elements that the matrix
* will contain. This is just an optional hint. If the actual number of non-zero
* elements turns out to be greater than the provided value, the initialization
* of the matrix will be slower because reallocation of internal data may be
* needed
*/
void Discretizer::initialize(long nUnknowns, long nMaximumNZ)
{
    _initialize(nUnknowns, nMaximumNZ);
}

/*!
* Internal function to initialize the discretizer.
*
* \param nUnknowns is the number of unknowns
* \param nMaximumNZ is the maximum number of non-zero elements that the matrix
* will contain. This is just an optional hint. If the actual number of non-zero
* elements turns out to be greater than the provided value, the initialization
* of the matrix will be slower because reallocation of internal data may be
* needed
*/
void Discretizer::_initialize(long nUnknowns, long nMaximumNZ)
{
    m_matrix.initialize(nUnknowns, nUnknowns, nMaximumNZ);
    m_constants.resize(nUnknowns);
}

/*!
 * Clear the discretizer
 */
void Discretizer::clear(bool release)
{
    SystemSolver::clear();

    m_matrix.clear(release);
    if (release) {
        std::vector<double>().swap(m_constants);
    } else {
        m_constants.clear();
    }
}

/**
* Add a stencil.
*
* \param stencil is the stencil that wil be added
*/
void Discretizer::importStencil(const StencilScalar &stencil)
{
    m_matrix.addRow(stencil.getPattern().getItemCount(), stencil.getPattern().data(), stencil.getWeights().data());
    m_constants.push_back(stencil.getConstant());
}

/*!
 * Solve the system
 */
void Discretizer::solve()
{
    // Initialzie the system
    if (isInitialized()) {
        initialize(m_matrix.getRowCount(), m_matrix.getNZCount());
    }

    // Subtract constant terms to the RHS
    long nUnknowns = m_matrix.getRowCount();
    double *raw_rhs = getRHSRawPtr();
    for (int i = 0; i < nUnknowns; ++i) {
        raw_rhs[i] -= m_constants[i];
    }
    restoreRHSRawPtr(raw_rhs);

    // Solve the system
    SystemSolver::solve();
}

/*!
 * Solve the system
 *
 * \param rhs is the right-hand-side of the system
 * \param solution in input should contain the initial solution, on output it
 * contains the solution of the linear system
 */
void Discretizer::solve(const std::vector<double> &rhs, std::vector<double> *solution)
{
    // Fills the vectors
    vectorsFill(rhs, solution);

    // Solve the system
    solve();

    // Export the solution
    vectorsExport(solution);
}

}
