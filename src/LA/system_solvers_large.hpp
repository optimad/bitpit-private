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

#ifndef __BITPIT_SYSTEM_SOLVERS_LARGE_HPP__
#define __BITPIT_SYSTEM_SOLVERS_LARGE_HPP__

#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <petscksp.h>

#include "system_matrix.hpp"

namespace bitpit {

struct KSPOptions {
    PetscInt restart;
    PetscInt levels;
    PetscInt overlap;
    PetscInt sublevels;
    PetscInt maxits;
    PetscScalar rtol;
    PetscScalar subrtol;
    bool nullspace;

    KSPOptions()
        : restart(50), levels(1), overlap(0), sublevels(4),
        maxits(10000), rtol(1.e-13), subrtol(1.e-13), nullspace(false)
    {
    }
};

struct KSPStatus {
    PetscErrorCode error;
    PetscInt its;
    KSPConvergedReason convergence;

    KSPStatus()
        : error(0), its(-1), convergence(KSP_DIVERGED_BREAKDOWN)
    {
    }
};

class SystemSolver {

public:
    enum PivotType {
        PIVOT_NONE,  // Natural
        PIVOT_ND,    // Nested Dissection
        PIVOT_1WD,   // One-way Dissection
        PIVOT_RCM,   // Reverse Cuthill-McKee
        PIVOT_MD     // Quotient Minimum Degree
    };

    static void addInitOption(std::string options);
    static void addInitOptions(const std::vector<std::string> &options);

#if BITPIT_ENABLE_MPI==1
    SystemSolver(MPI_Comm comm, bool debug = false);
#else
    SystemSolver(bool debug = false);
#endif
    ~SystemSolver();

    void clear();
    void initialize(const SparseMatrix &matrix, PivotType pivotType = PIVOT_NONE);
    bool isInitialized() const;

    void solve();
    void solve(const std::vector<double> &rhs, std::vector<double> *solution);

    void dump(const std::string &directory, const std::string &prefix = "") const;

    PivotType getPivotType();

    KSPOptions & getKSPOptions();
    const KSPOptions & getKSPOptions() const;
    const KSPStatus & getKSPStatus() const;

    double * getRHSRawPtr();
    const double * getRHSRawPtr() const;
    const double * getRHSRawReadPtr() const;
    void restoreRHSRawPtr(double *raw_rhs);
    void restoreRHSRawReadPtr(const double *raw_rhs) const;

    double * getSolutionRawPtr();
    const double * getSolutionRawPtr() const;
    const double * getSolutionRawReadPtr() const;
    void restoreSolutionRawPtr(double *raw_solution);
    void restoreSolutionRawReadPtr(const double *raw_solution) const;

protected:
    void matrixInit(const SparseMatrix &matrix);
    void matrixFill(const SparseMatrix &matrix);
    void matrixReorder();

#if BITPIT_ENABLE_MPI == 1
    void vectorsInit(const std::vector<long> &ghosts);
#else
    void vectorsInit();
#endif
    void vectorsReorder(PetscBool inv);
    void vectorsFill(const std::vector<double> &rhs, std::vector<double> *solution);
    void vectorsExport(std::vector<double> *solution);

private:
    static int m_nInstances;
    static std::vector<std::string> m_options;

    bool m_initialized;
    PivotType m_pivotType;

#if BITPIT_ENABLE_MPI==1
    MPI_Comm m_communicator;
    long m_rowGlobalIdOffset;
#endif

    Mat m_A;
    Vec m_rhs;
    Vec m_solution;

    IS m_rpivot;
    IS m_cpivot;

    KSP m_KSP;
    KSPOptions m_KSPOptions;
    KSPStatus m_KSPStatus;

    void pivotInit(PivotType pivotType);

    void KSPInit();

};

}

#endif
