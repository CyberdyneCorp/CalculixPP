#include "calculixpp/numerics/linear_static.hpp"

#include <cstdint>

#include "numpp/core/creation.hpp"
#include "numpp/core/dtype.hpp"
#include "scipp/sparse/sparse.hpp"

namespace cxpp::numerics {

std::vector<Real> solve_reduced(const fem::LinearSystem& sys, SolverKind kind) {
  const std::int64_t n = sys.n_free;
  const std::int64_t nnz = static_cast<std::int64_t>(sys.vals.size());

  scipp::sparse::CooMatrix coo;
  coo.rows = n;
  coo.cols = n;
  coo.data = numpp::zeros({nnz}, numpp::kFloat64);
  coo.row = numpp::zeros({nnz}, numpp::kInt64);
  coo.col = numpp::zeros({nnz}, numpp::kInt64);
  for (std::int64_t i = 0; i < nnz; ++i) {
    coo.data.set_item<double>({i}, sys.vals[static_cast<std::size_t>(i)]);
    coo.row.set_item<std::int64_t>({i}, sys.rows[static_cast<std::size_t>(i)]);
    coo.col.set_item<std::int64_t>({i}, sys.cols[static_cast<std::size_t>(i)]);
  }

  const scipp::sparse::CsrMatrix A = scipp::sparse::CsrMatrix::from_coo(coo);

  numpp::ndarray b = numpp::zeros({n}, numpp::kFloat64);
  for (std::int64_t i = 0; i < n; ++i)
    b.set_item<double>({i}, sys.rhs[static_cast<std::size_t>(i)]);

  const numpp::ndarray x = (kind == SolverKind::CG)
                               ? scipp::sparse::cg(A, b)
                               : scipp::sparse::spsolve(A, b);

  std::vector<Real> u(static_cast<std::size_t>(n));
  for (std::int64_t i = 0; i < n; ++i)
    u[static_cast<std::size_t>(i)] = x.item<double>({i});
  return u;
}

LinearStaticResult solve_linear_static(const Model& model, SolverKind kind) {
  const fem::LinearSystem sys = fem::assemble_linear_static(model);
  const std::vector<Real> uf = solve_reduced(sys, kind);

  LinearStaticResult res;
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
  return res;
}

}  // namespace cxpp::numerics
