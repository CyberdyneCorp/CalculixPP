#include "calculixpp/numerics/linear_static.hpp"

#include "calculixpp/compute/backend.hpp"
#include "calculixpp/fem/stress.hpp"

namespace cxpp::numerics {

SolverKind solver_kind(const Model& model) {
  return model.solver == RequestedSolver::CG ? SolverKind::CG : SolverKind::Direct;
}

std::vector<Real> solve_reduced(const fem::LinearSystem& sys, SolverKind kind) {
  // Route the sparse SPD solve through the compute backend (CPU reference backend
  // wraps the SciPP path). Behavior-preserving: same COO triplets, same solver.
  const compute::ComputeBackend& backend = compute::select_backend();
  return backend.solve_sparse(sys.rows, sys.cols, sys.vals, sys.n_free, sys.rhs,
                              kind);
}

StaticFields solve_linear_static(const Model& model, SolverKind kind) {
  const fem::LinearSystem sys = fem::assemble_linear_static(model);
  const std::vector<Real> uf = solve_reduced(sys, kind);

  StaticFields res;
  const std::size_t n_nodes = model.mesh.num_nodes();
  res.displacement.assign(n_nodes, Vec3{0, 0, 0});
  for (std::size_t ni = 0; ni < n_nodes; ++ni) {
    for (int c = 0; c < kDofsPerNode; ++c) {
      const std::size_t g = ni * kDofsPerNode + static_cast<std::size_t>(c);
      const Index eq = sys.dof_eq[g];
      res.displacement[ni][static_cast<std::size_t>(c)] =
          (eq >= 0) ? uf[static_cast<std::size_t>(eq)] : sys.prescribed[g];
    }
  }
  fem::recover_fields(model, res);  // fills stress + reaction
  return res;
}

}  // namespace cxpp::numerics
