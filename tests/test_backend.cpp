// ComputeBackend + solver-selection tests.
//  - The CPU backend's solve_sparse reproduces the pre-backend solution on the
//    single-tet model (cross-checked against an independent dense solve).
//  - SOLVER= on *STATIC maps onto RequestedSolver; an unknown name raises.
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "calculixpp/compute/backend.hpp"
#include "calculixpp/fem/assembly.hpp"
#include "calculixpp/io/inp_parser.hpp"
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
    K[static_cast<std::size_t>(sys.rows[t])][static_cast<std::size_t>(sys.cols[t])] +=
        sys.vals[t];
  return dense_solve(K, sys.rhs);
}

void test_backend_matches_reference() {
  const Model m = make_single_tet();
  const fem::LinearSystem sys = fem::assemble_linear_static(m);
  const std::vector<Real> uref = reference_free_solution(sys);

  // select_backend() must yield the CPU backend regardless of the request; an
  // unimplemented backend falls back to CPU (never an error).
  const compute::ComputeBackend& cpu = compute::select_backend();
  CX_CHECK(cpu.kind() == compute::BackendKind::CPU);
  CX_CHECK(compute::select_backend(compute::BackendKind::CUDA).kind() ==
           compute::BackendKind::CPU);

  const std::vector<Real> u = cpu.solve_sparse(
      sys.rows, sys.cols, sys.vals, sys.n_free, sys.rhs, compute::SolverKind::Direct);
  CX_CHECK(u.size() == uref.size());
  for (std::size_t i = 0; i < u.size(); ++i)
    CX_NEAR(u[i], uref[i], 1e-8 * (1.0 + std::fabs(uref[i])));
}

// Minimal deck body with a swappable *STATIC card, for SOLVER= parsing.
std::string deck_with_static(const std::string& static_card) {
  return R"(*NODE
1, 0., 0., 0.
2, 1., 0., 0.
3, 0., 1., 0.
4, 0., 0., 1.
*ELEMENT, TYPE=C3D4, ELSET=EALL
100, 1, 2, 3, 4
*BOUNDARY
1, 1, 3
2, 1, 3
3, 1, 3
*MATERIAL, NAME=EL
*ELASTIC
210000., 0.3
*SOLID SECTION, ELSET=EALL, MATERIAL=EL
*STEP
)" + static_card + R"(
*CLOAD
4, 3, -1000.
*END STEP
)";
}

void test_solver_selection() {
  // Unspecified -> Direct.
  CX_CHECK(io::parse_inp(deck_with_static("*STATIC")).solver ==
           RequestedSolver::Direct);
  // Direct family.
  CX_CHECK(io::parse_inp(deck_with_static("*STATIC, SOLVER=SPOOLES")).solver ==
           RequestedSolver::Direct);
  CX_CHECK(io::parse_inp(deck_with_static("*STATIC, SOLVER=PARDISO")).solver ==
           RequestedSolver::Direct);
  CX_CHECK(io::parse_inp(deck_with_static("*STATIC, SOLVER=PASTIX")).solver ==
           RequestedSolver::Direct);
  // Iterative family -> CG.
  CX_CHECK(io::parse_inp(deck_with_static("*STATIC, SOLVER=ITERATIVE SCALING"))
               .solver == RequestedSolver::CG);
  CX_CHECK(io::parse_inp(deck_with_static("*STATIC, SOLVER=CG")).solver ==
           RequestedSolver::CG);

  // solver_kind maps the model request onto the numeric SolverKind.
  CX_CHECK(numerics::solver_kind(io::parse_inp(deck_with_static("*STATIC"))) ==
           numerics::SolverKind::Direct);
  CX_CHECK(numerics::solver_kind(
               io::parse_inp(deck_with_static("*STATIC, SOLVER=CG"))) ==
           numerics::SolverKind::CG);
}

void test_unknown_solver_raises() {
  bool threw = false;
  std::string msg;
  try {
    io::parse_inp(deck_with_static("*STATIC, SOLVER=BOGUS"));
  } catch (const io::ParseError& e) {
    threw = true;
    msg = e.what();
  }
  CX_CHECK(threw);
  CX_CHECK(msg.find("BOGUS") != std::string::npos);  // names the requested solver
}

void test_backend_introspection() {
  // Names round-trip through backend_kind (case-insensitive).
  CX_CHECK(compute::backend_name(compute::BackendKind::CPU) == "cpu");
  CX_CHECK(compute::backend_name(compute::BackendKind::CUDA) == "cuda");
  CX_CHECK(compute::backend_kind("CPU") == compute::BackendKind::CPU);
  CX_CHECK(compute::backend_kind("Metal") == compute::BackendKind::Metal);

  bool threw = false;
  try {
    compute::backend_kind("bogus");
  } catch (const std::invalid_argument&) {
    threw = true;
  }
  CX_CHECK(threw);

  // Only CPU is implemented in Phase 1.
  CX_CHECK(compute::backend_available(compute::BackendKind::CPU));
  CX_CHECK(!compute::backend_available(compute::BackendKind::CUDA));
  const std::vector<compute::BackendKind> avail = compute::available_backends();
  CX_CHECK(avail.size() == 1);
  CX_CHECK(avail.front() == compute::BackendKind::CPU);
}

}  // namespace

int main() {
  test_backend_matches_reference();
  test_solver_selection();
  test_unknown_solver_raises();
  test_backend_introspection();
  if (cxtest::g_failures == 0) std::printf("test_backend: OK\n");
  CX_MAIN_RETURN();
}
