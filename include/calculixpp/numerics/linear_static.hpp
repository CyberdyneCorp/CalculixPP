#pragma once
#include <vector>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/types.hpp"
#include "calculixpp/fem/assembly.hpp"

// Sparse linear-static solve. Wraps SciPP's sparse module (scipp::sparse): the
// assembled COO triplets become a CsrMatrix and are solved with spsolve (direct)
// or cg (SPD iterative). NumPP has no sparse module — sparse lives in SciPP.
// (spec: linear-algebra-and-solvers, static-analysis.)
namespace cxpp::numerics {

enum class SolverKind {
  Direct,  // scipp::sparse::spsolve (LU-based direct)
  CG,      // scipp::sparse::cg (SPD iterative)
};

// Solve the reduced free/free system; returns the free-DOF solution (size n_free).
std::vector<Real> solve_reduced(const fem::LinearSystem& sys,
                                SolverKind kind = SolverKind::Direct);

struct LinearStaticResult {
  // Nodal displacement, aligned with mesh node indices.
  std::vector<Vec3> displacement;
};

// Assemble and solve a linear-static step end to end.
LinearStaticResult solve_linear_static(const Model& model,
                                       SolverKind kind = SolverKind::Direct);

}  // namespace cxpp::numerics
