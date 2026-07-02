// Assembly + reduced-solve test on a single C3D4 tet: fix three nodes, load the
// fourth, assemble the 3x3 free system, solve it densely, and verify the residual.
#include <array>
#include <vector>

#include "calculixpp/fem/assembly.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

// Dense Gaussian elimination for the small reduced system (test-only).
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
  // Fully fix nodes 1,2,3.
  for (Index nd : {1, 2, 3})
    for (int c = 1; c <= 3; ++c) m.spcs.push_back(Spc{nd, c, 0.0});
  // Load node 4 in -z.
  m.cloads.push_back(Cload{4, 3, -1000.0});
  return m;
}

}  // namespace

int main() {
  const Model m = make_single_tet();
  const fem::LinearSystem sys = fem::assemble_linear_static(m);

  CX_CHECK(sys.n_free == 3);  // only node 4's 3 DOFs are free
  CX_CHECK(sys.rhs.size() == 3);

  // Rebuild dense K_ff from COO, check symmetry + positive diagonal.
  const std::size_t n = static_cast<std::size_t>(sys.n_free);
  std::vector<std::vector<Real>> K(n, std::vector<Real>(n, 0.0));
  for (std::size_t t = 0; t < sys.vals.size(); ++t)
    K[static_cast<std::size_t>(sys.rows[t])][static_cast<std::size_t>(sys.cols[t])] += sys.vals[t];
  for (std::size_t i = 0; i < n; ++i) {
    CX_CHECK(K[i][i] > 0.0);
    for (std::size_t j = 0; j < n; ++j) CX_NEAR(K[i][j], K[j][i], 1e-6);
  }

  // Solve and verify residual K u = f.
  const std::vector<Real> u = dense_solve(K, sys.rhs);
  for (std::size_t i = 0; i < n; ++i) {
    Real r = 0.0;
    for (std::size_t j = 0; j < n; ++j) r += K[i][j] * u[j];
    CX_NEAR(r, sys.rhs[i], 1e-6);
  }
  // The loaded node should displace in -z (u_z < 0) and be non-trivial.
  CX_CHECK(u[2] < 0.0);

  if (cxtest::g_failures == 0) std::printf("test_assembly: OK\n");
  CX_MAIN_RETURN();
}
