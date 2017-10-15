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

#ifndef __BITPIT_DISCRETIZER_HPP__
#define __BITPIT_DISCRETIZER_HPP__

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "bitpit_LA.hpp"

#include "discretizer.hpp"
#include "stencil.hpp"

namespace bitpit {

class Discretizer : protected SystemSolver {

public:
#if BITPIT_ENABLE_MPI==1
    Discretizer(MPI_Comm comm, bool debug = false);
    Discretizer(MPI_Comm comm, long nUnknowns, long nMaximumNZ = 0, bool debug = false);
#else
    Discretizer(bool debug = false);
    Discretizer(long nUnknowns, long nMaximumNZ = 0, bool debug = false);
#endif

    void initialize(long nUnknowns, long nMaximumNZ);
    void clear(bool release = false);

    void importStencil(const StencilScalar &stencil);

    void solve();
    void solve(const std::vector<double> &rhs, std::vector<double> *solution);

    using SystemSolver::getRHSRawPtr;
    using SystemSolver::getSolutionRawPtr;

protected:
    SparseMatrix m_matrix;
    std::vector<double> m_constants;

    void _initialize(long nUnknowns, long nMaximumNZ);

};

}

#endif
