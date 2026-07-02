// End-to-end sparse solve via SciPP, cross-checked against an independent dense
// solve on the single-tet model (fix nodes 1-3, load node 4 in -z).
#include <array>
#include <cmath>
#include <vector>

#include "calculixpp/fem/assembly.hpp"
#include "calculixpp/numerics/linear_static.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

std::vector<Real> dense_solve(std::vector<std::vector<Real>> A, std::vector<Real> b) {
  const std::size_t n = b.size();
  for (std::size_t c = 0; c < n; ++c) {
    std::size_t piv = c;
    for (std::size_t r = c + 1; r < n; ++r)
      if (std::fabs(A[r][c]) > std::fabs(A[piv][c])) piv = r;
    std::swap(A[c], A[piv]);
    std::swap(b[c], b[piv]);
    for (std::size_t r = c + 1; r < n; ++r) {
      const Real f = A[r][c] / A[c][c];
      for (std::size_t k = c; k < n; ++k) A[r][k] -= f * A[c][k];
      b[r] -= f * b[c];
    }
  }
  std::vector<Real> x(n, 0.0);
  for (std::size_t i = n; i-- > 0;) {
    Real s = b[i];
    for (std::size_t k = i + 1; k < n; ++k) s -= A[i][k] * x[k];
    x[i] = s / A[i][i];
  }
  return x;
}

Model make_single_tet() {
  Model m;
  m.mesh.add_node(1, {0, 0, 0});
  m.mesh.add_node(2, {1, 0, 0});
  m.mesh.add_node(3, {0, 1, 0});
  m.mesh.add_node(4, {0, 0, 1});
  m.mesh.add_element(100, ElementType::C3D4, {1, 2, 3, 4});
  m.mesh.add_elset("EALL", {100});
  m.materials["EL"] = Material{"EL", ElasticIso{210000.0, 0.3}, std::nullopt};
  m.sections.push_back(SolidSection{"EALL", "EL"});
  for (Index nd : {1, 2, 3})
    for (int c = 1; c <= 3; ++c) m.spcs.push_back(Spc{nd, c, 0.0});
  m.cloads.push_back(Cload{4, 3, -1000.0});
  return m;
}

std::vector<Real> reference_free_solution(const fem::LinearSystem& sys) {
  const std::size_t n = static_cast<std::size_t>(sys.n_free);
  std::vector<std::vector<Real>> K(n, std::vector<Real>(n, 0.0));
  for (std::size_t t = 0; t < sys.vals.size(); ++t)
    K[static_cast<std::size_t>(sys.rows[t])][static_cast<std::size_t>(sys.cols[t])] += sys.vals[t];
  return dense_solve(K, sys.rhs);
}

}  // namespace

int main() {
  const Model m = make_single_tet();
  const fem::LinearSystem sys = fem::assemble_linear_static(m);
  const std::vector<Real> uref = reference_free_solution(sys);

  for (numerics::SolverKind kind :
       {numerics::SolverKind::Direct, numerics::SolverKind::CG}) {
    const std::vector<Real> uf = numerics::solve_reduced(sys, kind);
    CX_CHECK(uf.size() == uref.size());
    for (std::size_t i = 0; i < uf.size(); ++i)
      CX_NEAR(uf[i], uref[i], 1e-8 * (1.0 + std::fabs(uref[i])));
  }

  // Full pipeline: node 4 (index 3) displaces in -z and matches the reference.
  const numerics::LinearStaticResult res = numerics::solve_linear_static(m);
  CX_CHECK(res.displacement.size() == 4);
  CX_CHECK(res.displacement[3][2] < 0.0);
  CX_NEAR(res.displacement[3][0], uref[0], 1e-8 * (1.0 + std::fabs(uref[0])));
  CX_NEAR(res.displacement[3][1], uref[1], 1e-8 * (1.0 + std::fabs(uref[1])));
  CX_NEAR(res.displacement[3][2], uref[2], 1e-8 * (1.0 + std::fabs(uref[2])));

  if (cxtest::g_failures == 0) std::printf("test_solve: OK\n");
  CX_MAIN_RETURN();
}
