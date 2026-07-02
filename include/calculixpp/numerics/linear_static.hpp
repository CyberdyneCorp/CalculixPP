#pragma once
#include <vector>

#include "calculixpp/compute/backend.hpp"
#include "calculixpp/core/model.hpp"
#include "calculixpp/core/results.hpp"
#include "calculixpp/core/types.hpp"
#include "calculixpp/fem/assembly.hpp"

// Sparse linear-static solve. Routed through the ComputeBackend abstraction; the
// CPU reference backend wraps SciPP's sparse module (scipp::sparse): the assembled
// COO triplets become a CsrMatrix and are solved with spsolve (direct) or cg (SPD
// iterative). NumPP has no sparse module — sparse lives in SciPP.
// (spec: linear-algebra-and-solvers, static-analysis, compute-backend.)
namespace cxpp::numerics {

// PERFORMANCE CAVEAT (SciPP#10): scipp::sparse::spsolve currently densifies the
// matrix (O(N^2) mem / O(N^3) time) — correct but caps mesh size at ~hundreds of
// DOF; scipp::sparse::cg is unpreconditioned and may NOT converge on stiff FE
// systems (silently wrong). Revisit defaults once SciPP ships a sparse Cholesky/
// LDLT + preconditioned CG. See https://github.com/CyberdyneCorp/SciPP/issues/10
//
// SolverKind lives in the compute layer (compute/backend.hpp) so the backend
// interface and numerics share one enum; re-exported here for existing callers.
using SolverKind = compute::SolverKind;

// Map the model's requested solver (from SOLVER= on *STATIC) to a SolverKind.
SolverKind solver_kind(const Model& model);

// Solve the reduced free/free system; returns the free-DOF solution (size n_free).
std::vector<Real> solve_reduced(const fem::LinearSystem& sys,
                                SolverKind kind = SolverKind::Direct);

// Assemble and solve a linear-static step end to end, returning nodal
// displacement, averaged stress, and reaction forces.
StaticFields solve_linear_static(const Model& model,
                                 SolverKind kind = SolverKind::Direct);

}  // namespace cxpp::numerics
