#pragma once
#include <vector>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/results.hpp"
#include "calculixpp/core/types.hpp"
#include "calculixpp/fem/assembly.hpp"

// Sparse linear-static solve. Wraps SciPP's sparse module (scipp::sparse): the
// assembled COO triplets become a CsrMatrix and are solved with spsolve (direct)
// or cg (SPD iterative). NumPP has no sparse module — sparse lives in SciPP.
// (spec: linear-algebra-and-solvers, static-analysis.)
namespace cxpp::numerics {

// PERFORMANCE CAVEAT (SciPP#10): scipp::sparse::spsolve currently densifies the
// matrix (O(N^2) mem / O(N^3) time) — correct but caps mesh size at ~hundreds of
// DOF; scipp::sparse::cg is unpreconditioned and may NOT converge on stiff FE
// systems (silently wrong). Revisit defaults once SciPP ships a sparse Cholesky/
// LDLT + preconditioned CG. See https://github.com/CyberdyneCorp/SciPP/issues/10
enum class SolverKind {
  Direct,  // scipp::sparse::spsolve (dense LU today — see caveat)
  CG,      // scipp::sparse::cg (unpreconditioned — see caveat)
};

// Solve the reduced free/free system; returns the free-DOF solution (size n_free).
std::vector<Real> solve_reduced(const fem::LinearSystem& sys,
                                SolverKind kind = SolverKind::Direct);

// Assemble and solve a linear-static step end to end, returning nodal
// displacement, averaged stress, and reaction forces.
StaticFields solve_linear_static(const Model& model,
                                 SolverKind kind = SolverKind::Direct);

}  // namespace cxpp::numerics
