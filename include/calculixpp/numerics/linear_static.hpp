#pragma once
#include <optional>
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

// Solver notes (SciPP#10, resolved in SciPP v1.2.0): scipp::sparse::spsolve now
// factors genuinely sparsely (SPD -> Cholesky, else LU, with RCM fill-reducing
// ordering) — an ~280x speedup at 8k DOF and no longer memory-bound. The CG path
// (compute/cpu_backend.cpp) uses IC0-preconditioned cg_report and raises on
// non-convergence. Direct remains the default.
//
// SolverKind lives in the compute layer (compute/backend.hpp) so the backend
// interface and numerics share one enum; re-exported here for existing callers.
using SolverKind = compute::SolverKind;

// DOF count above which the Auto policy prefers IC0-CG over sparse-direct.
// Heuristic: sparse Cholesky is fastest and exact for small/medium systems, but
// its factor fill-in grows super-linearly for large 3D meshes and can exceed
// device memory (esp. mobile); above the threshold, IC0-CG stays O(nnz) in
// memory. Tunable per platform. (Direct is still faster at ~10^4 DOF — the switch
// trades speed for memory scalability at large N.)
inline constexpr Index kAutoDirectMaxDof = 50000;

// Resolve a requested solver + system size to a concrete SolverKind. Auto ->
// Direct for n_free <= direct_max_dof, else CG; Direct/CG pass through.
SolverKind resolve_solver_kind(RequestedSolver req, Index n_free,
                               Index direct_max_dof = kAutoDirectMaxDof);

// Solve the reduced free/free system; returns the free-DOF solution (size n_free).
std::vector<Real> solve_reduced(const fem::LinearSystem& sys,
                                SolverKind kind = SolverKind::Direct);

// Assemble and solve a linear-static step end to end, returning nodal
// displacement, averaged stress, and reaction forces. When `forced` is empty the
// solver is chosen from the model's SOLVER= request via resolve_solver_kind
// (Auto uses the system size); pass a value to override.
StaticFields solve_linear_static(const Model& model,
                                 std::optional<SolverKind> forced = std::nullopt);

}  // namespace cxpp::numerics
