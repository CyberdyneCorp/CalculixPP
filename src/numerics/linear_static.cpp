#include "calculixpp/numerics/linear_static.hpp"

#include "calculixpp/compute/backend.hpp"
#include "calculixpp/fem/stress.hpp"

namespace cxpp::numerics {

SolverKind resolve_solver_kind(RequestedSolver req, Index n_free,
                               Index direct_max_dof) {
  switch (req) {
    case RequestedSolver::Direct:
      return SolverKind::Direct;
    case RequestedSolver::CG:
      return SolverKind::CG;
    case RequestedSolver::Auto:
    default:
      return n_free > direct_max_dof ? SolverKind::CG : SolverKind::Direct;
  }
}

std::vector<Real> solve_reduced(const fem::LinearSystem& sys, SolverKind kind) {
  // Route the sparse SPD solve through the compute backend (CPU reference backend
  // wraps the SciPP path). Behavior-preserving: same COO triplets, same solver.
  const compute::ComputeBackend& backend = compute::select_backend();
  return backend.solve_sparse(sys.rows, sys.cols, sys.vals, sys.n_free, sys.rhs,
                              kind);
}

StaticFields solve_linear_static(const Model& model,
                                 std::optional<SolverKind> forced) {
  const fem::LinearSystem sys = fem::assemble_linear_static(model);
  const SolverKind kind =
      forced ? *forced : resolve_solver_kind(model.solver, sys.n_free);
  const std::vector<Real> uf = solve_reduced(sys, kind);

  StaticFields res;
  const std::size_t n_nodes = model.mesh.num_nodes();
  res.displacement.assign(n_nodes, Vec3{0, 0, 0});
  // Expand the free solution to full nodal displacement through the constraint
  // transform: free DOFs take uf[eq], SPC DOFs take their prescribed value (full
  // magnitude here), MPC slaves are reconstructed from their masters.
  for (std::size_t ni = 0; ni < n_nodes; ++ni)
    for (int c = 0; c < kDofsPerNode; ++c) {
      const std::size_t g = ni * kDofsPerNode + static_cast<std::size_t>(c);
      res.displacement[ni][static_cast<std::size_t>(c)] =
          sys.transform.displacement(g, uf, sys.prescribed);
    }
  fem::recover_fields(model, res);  // fills stress + reaction
  return res;
}

}  // namespace cxpp::numerics
