// Damped complex modal eigensolver (*COMPLEX FREQUENCY, proportional damping — option B).
//   * ANALYTICAL Rayleigh-damped SDOF: the reduced complex eigenvalue equals the exact
//     closed form λ = -ζω ± iω√(1-ζ²) with ζ = (α/ω + β·ω)/2 (task 7.1);
//   * undamped limit α=β=0 reproduces the real *FREQUENCY ω with ζ=0 (task 7.3);
//   * a small Rayleigh-damped 2-DOF chain matches the per-mode closed form (task 7.2);
//   * the modal-reduced path agrees with the dense full 2n state-space oracle (task 7.4).
// The validation is analytical (exact for proportional damping) — deliberately NOT a
// CalculiX *COMPLEX FREQUENCY reference comparison, which solves the gyroscopic problem.
#include <cmath>
#include <complex>
#include <vector>

#include "calculixpp/fem/assembly.hpp"
#include "calculixpp/numerics/eigensolution.hpp"
#include "calculixpp/numerics/modal_dynamics.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

// Exact under-damped SDOF eigenvalue for the Im>=0 branch: λ = -ζω + iω√(1-ζ²).
std::complex<Real> closed_form(Real omega, Real zeta) {
  return {-zeta * omega, omega * std::sqrt(std::max(1.0 - zeta * zeta, 0.0))};
}

// A single grounded spring-mass on x: ω = sqrt(k/m), one free DOF.
Model make_sdof(Real k, Real mass) {
  Model m;
  m.mesh.add_node(1, {0, 0, 0});
  Spring s;
  s.kind = Spring::Kind::Grounded;
  s.node1 = 1;
  s.dof1 = 1;
  s.stiffness = k;
  m.springs = {s};
  m.point_masses = {PointMass{1, mass}};
  for (int c : {2, 3}) m.spcs.push_back(Spc{1, c, 0.0});  // free only on x
  m.procedure = Procedure::ComplexFrequency;
  return m;
}

// 2-DOF grounded chain on x: [[2k,-k],[-k,k]] over diag(m,m), analytical ω² known.
Model make_chain(Real k, Real mass) {
  Model m;
  m.mesh.add_node(1, {0, 0, 0});
  m.mesh.add_node(2, {1, 0, 0});
  Spring s1;
  s1.kind = Spring::Kind::Grounded;
  s1.node1 = 1;
  s1.dof1 = 1;
  s1.stiffness = k;
  Spring s2;
  s2.kind = Spring::Kind::Dof;
  s2.node1 = 1;
  s2.dof1 = 1;
  s2.node2 = 2;
  s2.dof2 = 1;
  s2.stiffness = k;
  m.springs = {s1, s2};
  m.point_masses = {PointMass{1, mass}, PointMass{2, mass}};
  for (Index n : {Index{1}, Index{2}})
    for (int c : {2, 3}) m.spcs.push_back(Spc{n, c, 0.0});
  m.procedure = Procedure::ComplexFrequency;
  return m;
}

// (7.1) Rayleigh-damped SDOF: the single complex eigenvalue equals the closed form to
// ~1e-8 relative, and the derived scalars (ω_d, ω_n, ζ, f_d) are consistent.
void test_sdof_closed_form() {
  const Real k = 400.0, mass = 4.0;  // ω = sqrt(100) = 10
  Model m = make_sdof(k, mass);
  m.rayleigh = Rayleigh{0.5, 0.002, true};  // ζ = (0.5/10 + 0.002*10)/2 = 0.035

  const fem::LinearSystem K = fem::assemble_linear_static(m);
  const fem::LinearSystem M = fem::assemble_mass(m, false);
  const numerics::EigenBasis basis = numerics::extract_modes(K, M, 1);
  const Real omega = basis.modes[0].omega;
  CX_NEAR(omega, 10.0, 1e-9);

  numerics::Damping damp;
  damp.alpha = m.rayleigh.alpha;
  damp.beta = m.rayleigh.beta;
  const numerics::ComplexEigenBasis cx = numerics::extract_complex_modes(basis, damp, 1);
  CX_CHECK(cx.modes.size() == 1);

  const Real zeta = damp.ratio(0, omega);
  const std::complex<Real> exact = closed_form(omega, zeta);
  const numerics::ComplexMode& md = cx.modes[0];
  CX_NEAR(md.eigenvalue.real(), exact.real(), 1e-8 * std::abs(exact));
  CX_NEAR(md.eigenvalue.imag(), exact.imag(), 1e-8 * std::abs(exact));
  CX_NEAR(md.omega_n, omega, 1e-8 * omega);
  CX_NEAR(md.zeta, zeta, 1e-8);
  CX_NEAR(md.omega_d, omega * std::sqrt(1.0 - zeta * zeta), 1e-8 * omega);
  CX_NEAR(md.frequency, md.omega_d / (2.0 * M_PI), 1e-12);
}

// (7.3) Undamped limit α=β=0 reproduces the real *FREQUENCY ω with ζ=0.
void test_undamped_limit() {
  Model m = make_sdof(400.0, 4.0);  // ω = 10
  const fem::LinearSystem K = fem::assemble_linear_static(m);
  const fem::LinearSystem M = fem::assemble_mass(m, false);
  const numerics::EigenBasis basis = numerics::extract_modes(K, M, 1);
  numerics::Damping damp;  // α=β=0, no modal ratios -> undamped
  const numerics::ComplexEigenBasis cx = numerics::extract_complex_modes(basis, damp, 1);
  CX_CHECK(cx.modes.size() == 1);
  CX_NEAR(cx.modes[0].zeta, 0.0, 1e-10);
  CX_NEAR(cx.modes[0].omega_d, 10.0, 1e-8);
  CX_NEAR(cx.modes[0].decay_rate, 0.0, 1e-8);
}

// (7.2) Rayleigh-damped 2-DOF chain: each mode's complex eigenvalue matches its own
// closed form λ_k = -ζ_k ω_k + i ω_k √(1-ζ_k²).
void test_chain_closed_form() {
  const Real k = 100.0, mass = 2.0;  // r = k/m = 50
  Model m = make_chain(k, mass);
  m.rayleigh = Rayleigh{0.2, 0.001, true};

  const fem::LinearSystem K = fem::assemble_linear_static(m);
  const fem::LinearSystem M = fem::assemble_mass(m, false);
  const numerics::EigenBasis basis = numerics::extract_modes(K, M, 2);
  CX_CHECK(basis.modes.size() == 2);

  numerics::Damping damp;
  damp.alpha = m.rayleigh.alpha;
  damp.beta = m.rayleigh.beta;
  const numerics::ComplexEigenBasis cx = numerics::extract_complex_modes(basis, damp, 2);
  CX_CHECK(cx.modes.size() == 2);
  for (std::size_t j = 0; j < 2; ++j) {
    const Real omega = basis.modes[j].omega;
    const Real zeta = damp.ratio(j, omega);
    const std::complex<Real> exact = closed_form(omega, zeta);
    CX_NEAR(cx.modes[j].eigenvalue.real(), exact.real(), 1e-8 * std::abs(exact));
    CX_NEAR(cx.modes[j].eigenvalue.imag(), exact.imag(), 1e-8 * std::abs(exact));
  }
}

// (7.4) Cross-check: the modal-reduced complex modes agree with the dense full 2n
// state-space oracle [[0,I],[-M⁻¹K,-M⁻¹C]] for a small problem. Because Rayleigh damping
// keeps C_r diagonal, the reduction is EXACT, so the agreement is tight.
void test_dense_oracle_cross_check() {
  const Real k = 100.0, mass = 2.0;
  Model m = make_chain(k, mass);
  m.rayleigh = Rayleigh{0.3, 0.0015, true};

  const fem::LinearSystem K = fem::assemble_linear_static(m);
  const fem::LinearSystem M = fem::assemble_mass(m, false);
  const numerics::EigenBasis basis = numerics::extract_modes(K, M, 2);
  numerics::Damping damp;
  damp.alpha = m.rayleigh.alpha;
  damp.beta = m.rayleigh.beta;

  const numerics::ComplexEigenBasis reduced =
      numerics::extract_complex_modes(basis, damp, 2);
  const numerics::ComplexEigenBasis dense =
      numerics::extract_complex_modes_dense(K, M, damp, 2);
  CX_CHECK(reduced.modes.size() == dense.modes.size());
  for (std::size_t j = 0; j < reduced.modes.size(); ++j) {
    const Real scale = std::max(std::abs(dense.modes[j].eigenvalue), 1.0);
    CX_NEAR(reduced.modes[j].eigenvalue.real(), dense.modes[j].eigenvalue.real(),
            1e-9 * scale);
    CX_NEAR(reduced.modes[j].eigenvalue.imag(), dense.modes[j].eigenvalue.imag(),
            1e-9 * scale);
  }
}

}  // namespace

int main() {
  test_sdof_closed_form();
  test_undamped_limit();
  test_chain_closed_form();
  test_dense_oracle_cross_check();
  CX_MAIN_RETURN();
}
