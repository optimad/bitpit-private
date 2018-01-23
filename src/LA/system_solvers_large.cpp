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

#include "system_solvers_large.hpp"

namespace bitpit {

int SystemSolver::m_nInstances = 0;
std::vector<std::string> SystemSolver::m_options = std::vector<std::string>(1, "bitpit");

/*!
 * Set initialization option
 */
void SystemSolver::addInitOption(std::string option)
{
    if (m_nInstances != 0) {
        throw std::runtime_error("Initialization opions can be set only before initializing the solver.");
    }

    m_options.push_back(option);
}

/*!
 * Set initialization options
 */
void SystemSolver::addInitOptions(const std::vector<std::string> &options)
{
    if (m_nInstances != 0) {
        throw std::runtime_error("Initialization opions can be set only before initializing the solver.");
    }

    // The first option is the executable name and we set it to a dummy value.
    for (std::string option : options) {
        m_options.push_back(option);
    }
}

/*!
 * Default constuctor
 */
#if BITPIT_ENABLE_MPI==1
SystemSolver::SystemSolver(MPI_Comm communicator, bool debug)
#else
SystemSolver::SystemSolver(bool debug)
#endif
    : m_initialized(false), m_pivotType(PIVOT_NONE)
{
    // Add debug options
    if (debug) {
        addInitOption("-log_view");
        addInitOption("-ksp_monitor_true_residual");
        addInitOption("-ksp_converged_reason");
    }

    // Initialize Petsc
    if (m_nInstances == 0) {
        const char help[] = "None";

        int argc = m_options.size();
        char **argv = new char*[argc];
        for (int i = 0; i < argc; ++i) {
            argv[i] = (char*) m_options[i].c_str();
        }

        PetscInitialize(&argc, &argv, 0, help);

        delete[] argv;
    }

    // Increase the number of instances
    ++m_nInstances;

#if BITPIT_ENABLE_MPI==1
    // Create a communicator
    MPI_Comm_dup(communicator, &m_communicator);
#else
    // Set a fake MPI communicator
    m_communicator = PETSC_COMM_SELF;
#endif
}

/*!
 * Destructor
 */
SystemSolver::~SystemSolver()
{
    // Decrease the number of instances
    --m_nInstances;

#if BITPIT_ENABLE_MPI==1
    // Free the MPI communicator
    int finalizedCalled;
    MPI_Finalized(&finalizedCalled);
    if (!finalizedCalled) {
        MPI_Comm_free(&m_communicator);
    }
#endif

    // Finalize petsc
    if (m_nInstances == 0) {
        PetscFinalize();
    }
}

/*!
 * Clear the system
 */
void SystemSolver::clear()
{
    if (!m_initialized) {
        return;
    }

    MatDestroy(&m_A);
    VecDestroy(&m_rhs);
    VecDestroy(&m_solution);

    KSPDestroy(&m_KSP);

    if (getPivotType() != PIVOT_NONE) {
        ISDestroy(&m_rpivot);
        ISDestroy(&m_cpivot);
    }

#if BITPIT_ENABLE_MPI == 1
    m_rowGlobalIdOffset = 0;
#endif

    m_initialized = false;
}

/*!
 * Initialize the system.
 *
 * \param matrix is the matrix
 * \param pivotType is the type of pivoting that will be used
 */
void SystemSolver::initialize(const SparseMatrix &matrix, PivotType pivotType)
{
    // Check if the matrix is finalized
    if (!matrix.isFinalized()) {
        throw std::runtime_error("Unable to initialize the system. The matrix is not yet finalized..");
    }

    // Clear the system
    clear();

    // Initialize matrix
    matrixInit(matrix);
    matrixFill(matrix);

    // Initialize pivot
    pivotInit(pivotType);
    if (getPivotType() != PIVOT_NONE) {
        matrixReorder();
    }

    // Initialize RHS and solution vectors
#if BITPIT_ENABLE_MPI == 1
    vectorsInit(matrix.extractGhostGlobalCols());
#else
    vectorsInit();
#endif

    // Initialize Krylov solver
    KSPInit();

    // The system is now initialized
    m_initialized = true;
}

/*!
 * Check if the system is initialized.
 *
 * \return Returns true if the system is initialized, false otherwise.
 */
bool SystemSolver::isInitialized() const
{
    return m_initialized;
}

/*!
 * Solve the system
 */
void SystemSolver::solve()
{
    // Reorder the vectors
    if (getPivotType() != PIVOT_NONE) {
        vectorsReorder(PETSC_FALSE);
    }

    // Solve the system
    m_KSPStatus.error = KSPSolve(m_KSP, m_rhs, m_solution);

    // Set solver info
    if (m_KSPStatus.error == 0) {
        KSPGetIterationNumber(m_KSP, &m_KSPStatus.its);
        KSPGetConvergedReason(m_KSP, &m_KSPStatus.convergence);
    } else {
        m_KSPStatus.its         = -1;
        m_KSPStatus.convergence = KSP_DIVERGED_BREAKDOWN;
    }

    // Reorder the vectors
    if (getPivotType() != PIVOT_NONE) {
        vectorsReorder(PETSC_TRUE);
    }
}

/*!
 * Solve the system
 *
 * \param rhs is the right-hand-side of the system
 * \param solution in input should contain the initial solution, on output it
 * contains the solution of the linear system
 */
void SystemSolver::solve(const std::vector<double> &rhs, std::vector<double> *solution)
{
    // Fills the vectors
    vectorsFill(rhs, solution);

    // Solve the system
    solve();

    // Export the solution
    vectorsExport(solution);
}

/*!
 * Initializes the matrix.
 *
 * \param matrix is the matrix
 */
void SystemSolver::matrixInit(const SparseMatrix &matrix)
{
    long nRows = matrix.getRowCount();

    // Evaluate the offset for the numbering on this partition
    m_rowGlobalIdOffset = 0;

#if BITPIT_ENABLE_MPI == 1
    int nProcessors;
    MPI_Comm_size(m_communicator, &nProcessors);
    if (nProcessors > 1) {
        std::vector<long> nGlobalRows(nProcessors);
        MPI_Allgather(&nRows, 1, MPI_LONG, nGlobalRows.data(), 1, MPI_LONG, m_communicator);

        int rank;
        MPI_Comm_rank(m_communicator, &rank);
        for (int i = 0; i < rank; ++i) {
            m_rowGlobalIdOffset += nGlobalRows[i];
        }
    }
#endif

    // Evaluate the number of non-zero elements
    //
    // For each row we count the number of local non-zero elements (d_nnz) and
    // the number of non-zero elements that belong to other processors (o_nnz)
    std::vector<int> d_nnz(nRows);
    std::vector<int> o_nnz(nRows);

#if BITPIT_ENABLE_MPI == 1
    long firstRowGlobalId = m_rowGlobalIdOffset;
    long lastRowGlobalId  = firstRowGlobalId + nRows - 1;
#endif
    for (long row = 0; row < nRows; ++row) {
        ConstProxyVector<long> rowPattern = matrix.getRowPattern(row);
        int nRowNZ = rowPattern.size();

        d_nnz[row] = 0;
        o_nnz[row] = 0;
        for (int k = 0; k < nRowNZ; ++k) {
#if BITPIT_ENABLE_MPI == 1
            long columnGlobalId = rowPattern[k];
            if (columnGlobalId >= firstRowGlobalId && columnGlobalId <= lastRowGlobalId) {
                ++d_nnz[row];
            } else {
                ++o_nnz[row];
            }
#else
            ++d_nnz[row];
#endif
        }
    }

    // Create the matrix
#if BITPIT_ENABLE_MPI == 1
    MatCreateAIJ(m_communicator, nRows, nRows, PETSC_DETERMINE, PETSC_DETERMINE, 0, d_nnz.data(), 0, o_nnz.data(), &m_A);
#else
    MatCreateSeqAIJ(m_communicator, nRows, nRows, 0, d_nnz.data(), &m_A);
#endif
}

/*!
 * Fills the matrix.
 *
 * \param matrix is the matrix
 */
void SystemSolver::matrixFill(const SparseMatrix &matrix)
{
    const long nRows = matrix.getRowCount();
    const long maxRowNZ = matrix.getMaxRowNZCount();
    const long rowGlobalOffset = 0;

    // Create the matrix
    std::vector<PetscInt> rowNZGlobalIds(maxRowNZ);
    std::vector<PetscScalar> rowNZValues(maxRowNZ);

    for (long row = 0; row < nRows; ++row) {
        ConstProxyVector<long> rowPattern = matrix.getRowPattern(row);
        ConstProxyVector<double> rowValues = matrix.getRowValues(row);

        const int nRowNZ = rowPattern.size();
        const PetscInt globalRow = rowGlobalOffset + row;

        for (int k = 0; k < nRowNZ; ++k) {
            rowNZGlobalIds[k] = rowPattern[k];
            rowNZValues[k] = rowValues[k];
        }

        MatSetValues(m_A, 1, &globalRow, nRowNZ, rowNZGlobalIds.data(), rowNZValues.data(), INSERT_VALUES);
    }

    // Let petsc build the matrix
    MatAssemblyBegin(m_A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(m_A, MAT_FINAL_ASSEMBLY);
}

/*!
 * Reorder the matrix.
 */
void SystemSolver::matrixReorder()
{
    Mat B;
    MatPermute(m_A, m_rpivot, m_cpivot, &B);
    MatDestroy(&m_A);
    MatDuplicate(B, MAT_COPY_VALUES, &m_A);
    MatDestroy(&B);
}

/*!
 * Initialize rhs and solution vectors.
 */
#if BITPIT_ENABLE_MPI == 1
/*!
 * \param ghosts is the list of global ids that are ghosts for the local
 * processor
 */
void SystemSolver::vectorsInit(const std::vector<long> &ghosts)
#else
void SystemSolver::vectorsInit()
#endif
{
    PetscInt nRows;
    PetscInt nColumns;
    MatGetLocalSize(m_A, &nRows, &nColumns);

#if BITPIT_ENABLE_MPI == 1
    PetscInt nGlobalRows;
    PetscInt nGlobalColumns;
    MatGetSize(m_A, &nGlobalRows, &nGlobalColumns);

    int nGhosts = 0;
    std::vector<PetscInt> ghostGlobalIds(ghosts.size());
    for (auto &ghostGlobalId : ghosts) {
        ghostGlobalIds[nGhosts] = ghostGlobalId;
        ++nGhosts;
    }

    VecCreateGhost(m_communicator, nColumns, nGlobalColumns, nGhosts, ghostGlobalIds.data(), &m_solution);
    VecCreateGhost(m_communicator, nRows, nGlobalRows, nGhosts, ghostGlobalIds.data(), &m_rhs);
#else
    VecCreateSeq(m_communicator, nColumns, &m_solution);
    VecCreateSeq(m_communicator, nRows, &m_rhs);
#endif
}

/*!
 * Reorders rhs and solution vectors.
 */
void SystemSolver::vectorsReorder(PetscBool inv)
{
    VecPermute(m_solution, m_cpivot, inv);
    VecPermute(m_rhs, m_rpivot, inv);
}

/*!
 * Fills rhs and solution vectors.
 */
void SystemSolver::vectorsFill(const std::vector<double> &rhs, std::vector<double> *solution)
{
    // Import RHS
    int nRows;
    VecGetLocalSize(m_rhs, &nRows);

    PetscScalar *raw_rhs;
    VecGetArray(m_rhs, &raw_rhs);
    for (int i = 0; i < nRows; ++i) {
        raw_rhs[i] = rhs[i];
    }
    VecRestoreArray(m_rhs, &raw_rhs);

    // Import initial solution
    int nColumns;
    VecGetLocalSize(m_solution, &nColumns);

    PetscScalar *raw_solution;
    VecGetArray(m_solution, &raw_solution);
    for (int i = 0; i < nColumns; ++i) {
        raw_solution[i] = (*solution)[i];
    }
    VecRestoreArray(m_solution, &raw_solution);
}

/*!
 * Export the solution vector.
 */
void SystemSolver::vectorsExport(std::vector<double> *solution)
{
    int size;
    VecGetLocalSize(m_solution, &size);

    const PetscScalar *raw_solution;
    VecGetArrayRead(m_solution, &raw_solution);
    for (int i = 0; i < size; ++i) {
        (*solution)[i] = raw_solution[i];
    }
    VecRestoreArrayRead(m_solution, &raw_solution);
}

/*!
 * Get a raw pointer to the solution vector.
 *
 * \result A raw pointer to the solution vector.
 */
double * SystemSolver::getRHSRawPtr()
{
    PetscScalar *raw_rhs;
    VecGetArray(m_rhs, &raw_rhs);

    return raw_rhs;
}

/*!
 * Get a constant raw pointer to the solution vector.
 *
 * \result A constant raw pointer to the solution vector.
 */
const double * SystemSolver::getRHSRawPtr() const
{
    return getRHSRawReadPtr();
}

/*!
 * Get a constant raw pointer to the solution vector.
 *
 * \result A constant raw pointer to the solution vector.
 */
const double * SystemSolver::getRHSRawReadPtr() const
{
    const PetscScalar *raw_rhs;
    VecGetArrayRead(m_rhs, &raw_rhs);

    return raw_rhs;
}

/*!
 * Restores the solution vector after getRHSRawPtr() has been called.
 *
 * \param raw_solution is the location of pointer to array obtained from
 * getRHSRawPtr()
 */
void SystemSolver::restoreRHSRawPtr(double *raw_rhs)
{
    VecRestoreArray(m_rhs, &raw_rhs);
}

/*!
 * Restores the solution vector after getRHSRawReadPtr() has been called.
 *
 * \param raw_solution is the location of pointer to array obtained from
 * getRHSRawReadPtr()
 */
void SystemSolver::restoreRHSRawReadPtr(const double *raw_rhs) const
{
    VecRestoreArrayRead(m_rhs, &raw_rhs);
}

/*!
 * Get a raw pointer to the solution vector.
 *
 * \result A raw pointer to the solution vector.
 */
double * SystemSolver::getSolutionRawPtr()
{
    PetscScalar *raw_solution;
    VecGetArray(m_solution, &raw_solution);

    return raw_solution;
}

/*!
 * Get a constant raw pointer to the solution vector.
 *
 * \result A constant raw pointer to the solution vector.
 */
const double * SystemSolver::getSolutionRawPtr() const
{
    return getSolutionRawReadPtr();
}

/*!
 * Get a constant raw pointer to the solution vector.
 *
 * \result A constant raw pointer to the solution vector.
 */
const double * SystemSolver::getSolutionRawReadPtr() const
{
    const PetscScalar *raw_solution;
    VecGetArrayRead(m_solution, &raw_solution);

    return raw_solution;
}

/*!
 * Restores the solution vector after getSolutionRawPtr() has been called.
 *
 * \param raw_solution is the location of pointer to array obtained from
 * getSolutionRawPtr()
 */
void SystemSolver::restoreSolutionRawPtr(double *raw_solution)
{
    VecRestoreArray(m_solution, &raw_solution);
}

/*!
 * Restores the solution vector after getSolutionRawReadPtr() has been called.
 *
 * \param raw_solution is the location of pointer to array obtained from
 * getSolutionRawReadPtr()
 */
void SystemSolver::restoreSolutionRawReadPtr(const double *raw_solution) const
{
    VecRestoreArrayRead(m_solution, &raw_solution);
}

/*!
 * Dump the system to file
 */
void SystemSolver::dump(const std::string &directory, const std::string &prefix) const
{
    std::stringstream filePathStream;

    // Matrix
    PetscViewer matViewer;
    PetscViewerCreate(m_communicator, &matViewer);
    PetscViewerSetType(matViewer, PETSCVIEWERASCII);
    PetscViewerPushFormat(matViewer, PETSC_VIEWER_DEFAULT);

    filePathStream.str(std::string());
    filePathStream << directory << "/" << prefix << "A.txt";
    PetscViewerFileSetName(matViewer, filePathStream.str().c_str());
    MatView(m_A, matViewer);
    PetscViewerDestroy(&matViewer);

    // RHS
    PetscViewer rhsViewer;
    PetscViewerCreate(m_communicator, &rhsViewer);
    PetscViewerSetType(rhsViewer, PETSCVIEWERASCII);
    PetscViewerPushFormat(rhsViewer, PETSC_VIEWER_DEFAULT);

    filePathStream.str(std::string());
    filePathStream << directory << "/" << prefix << "rhs.txt";
    PetscViewerFileSetName(rhsViewer, filePathStream.str().c_str());
    VecView(m_rhs, rhsViewer);
    PetscViewerDestroy(&rhsViewer);
}

/*!
 * Initialize pivot
 *
 * \param pivotType is the type of pivoting that will be used
 */
void SystemSolver::pivotInit(PivotType pivotType)
{
    m_pivotType = pivotType;

    MatOrderingType petscPivotType;
    switch (m_pivotType) {

    case (PIVOT_ND):
        petscPivotType = MATORDERINGNATURAL;
        break;

    case (PIVOT_1WD):
        petscPivotType = MATORDERING1WD;
        break;

    case (PIVOT_RCM):
        petscPivotType = MATORDERINGRCM;
        break;

    case (PIVOT_MD):
        petscPivotType = MATORDERINGQMD;
        break;

    default:
        return;

    }

    MatGetOrdering(m_A, petscPivotType, &m_rpivot, &m_cpivot);
    ISSetPermutation(m_rpivot);
    ISSetPermutation(m_cpivot);
}

/*!
 * Get the pivot type.
 *
 * \result The pivot type.
 */
SystemSolver::PivotType SystemSolver::getPivotType()
{
    return m_pivotType;
}

/*!
 * Initialize the Krylov solver.
 */
void SystemSolver::KSPInit()
{
    int nProcessors;
#if BITPIT_ENABLE_MPI == 1
    MPI_Comm_size(m_communicator, &nProcessors);
#else
    nProcessors = 1;
#endif

    KSPCreate(m_communicator, &m_KSP);
    KSPSetOperators(m_KSP, m_A, m_A);

    if (m_KSPOptions.nullspace) {
        MatNullSpace nullspace;
        MatNullSpaceCreate(m_communicator, PETSC_TRUE, 0, NULL, &nullspace);
        MatSetNullSpace(m_A, nullspace);
        MatNullSpaceDestroy(&nullspace);
    }

    KSPSetType(m_KSP, KSPFGMRES);
    KSPGMRESSetRestart(m_KSP, m_KSPOptions.restart);
    KSPSetTolerances(m_KSP, m_KSPOptions.rtol, PETSC_DEFAULT, 1e10, m_KSPOptions.maxits);
    KSPSetInitialGuessNonzero(m_KSP, PETSC_TRUE);

    PC preconditioner;
    KSPGetPC(m_KSP, &preconditioner);
    if (nProcessors > 1) {
        PCSetType(preconditioner, PCASM);
        PCASMSetOverlap(preconditioner, m_KSPOptions.overlap);
    } else {
        PCSetType(preconditioner, PCILU);
        PCFactorSetLevels(preconditioner, m_KSPOptions.levels);
    }
    KSPSetFromOptions(m_KSP);
    KSPSetUp(m_KSP);

    // Set ASM sub block preconditioners
    if (nProcessors > 1) {
        KSP *subksp;
        PC subpc;
        PetscInt nlocal, first;
        PCASMGetSubKSP(preconditioner, &nlocal, &first, &subksp);
        for (PetscInt i = 0; i < nlocal; ++i) {
            KSPGetPC(subksp[i], &subpc);
            PCSetType(subpc, PCILU);
            PCFactorSetLevels(subpc, m_KSPOptions.sublevels);
            KSPSetTolerances(subksp[i], m_KSPOptions.subrtol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
        }
    }
}

/*!
 * Get a reference to the options associated to the Kryolov solver.
 *
 * \return A reference to the options associated to the Kryolov solver.
 */
KSPOptions & SystemSolver::getKSPOptions()
{
    return m_KSPOptions;
}

/*!
 * Get a constant reference to the options associated to the Kryolov solver.
 *
 * \return A constant reference to the options associated to the Kryolov solver.
 */
const KSPOptions & SystemSolver::getKSPOptions() const
{
    return m_KSPOptions;
}

/*!
 * Get a constant reference to the status of the Kryolov solver.
 *
 * \return A constant reference to the status of the Kryolov solver.
 */
const KSPStatus & SystemSolver::getKSPStatus() const
{
    return m_KSPStatus;
}

}
