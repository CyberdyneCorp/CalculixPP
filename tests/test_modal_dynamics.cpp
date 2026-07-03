// Modal superposition engine (Phase 4, tasks 1.5, 4.1, 4.2, 4.3).
// Validated ANALYTICALLY against closed-form SDOF / spring-mass results:
//   * free undamped modal response reproduces the exact period 2π/ω and amplitude
//     (the exact piecewise-linear integrator has zero algorithmic damping);
//   * a damped SDOF step response matches the closed-form under-damped solution;
//   * steady-state harmonic response: static deflection at Ω→0, resonant peak
//     F/k·Q (Q=1/2ζ), and correct half-power bandwidth Δω ≈ 2ζω;
//   * modal load projection p_k = φ_kᵀ f and base-motion effective load p_k = -Γ_k;
//   * Rayleigh damping maps to ζ_k = (α/ω_k + βω_k)/2.
#include <cmath>
#include <complex>
#include <vector>

#include "calculixpp/fem/assembly.hpp"
#include "calculixpp/numerics/eigensolution.hpp"
#include "calculixpp/numerics/modal_dynamics.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

// Build a single-DOF spring-mass system (mass m on node-1 x, grounded spring k), all
// other DOFs constrained, so K x = λ M x has the single mode ω = sqrt(k/m). Returns the
// eigenbasis (1 mode). This is the canonical analytical SDOF oscillator.
numerics::EigenBasis make_sdof(Real k, Real mass) {
  Model m;
  m.mesh.add_node(1, {0, 0, 0});
  Spring s;
  s.kind = Spring::Kind::Grounded;
  s.node1 = 1;
  s.dof1 = 1;
  s.stiffness = k;
  m.springs = {s};
  m.point_masses = {PointMass{1, mass}};
  for (int c : {2, 3}) m.spcs.push_back(Spc{1, c, 0.0});
  m.procedure = Procedure::Frequency;
  const fem::LinearSystem K = fem::assemble_linear_static(m);
  const fem::LinearSystem M = fem::assemble_mass(m, false);
  return numerics::extract_modes(K, M, 1);
}

// (4.1) Free undamped modal response: give the SDOF an initial modal displacement and
// integrate with zero load; the exact integrator must reproduce q(t) = q0 cos(ω t) with
// no amplitude decay over many periods, confirming the natural period 2π/ω.
void test_free_period_undamped() {
  const Real k = 400.0, mass = 1.0;  // ω = 20 rad/s
  numerics::EigenBasis basis = make_sdof(k, mass);
  const Real omega = basis.modes[0].omega;
  CX_NEAR(omega, std::sqrt(k / mass), 1e-9);

  numerics::Damping damp;  // undamped
  numerics::ModalSystem sys = numerics::project_modal_system(basis, damp);

  // Integrate the bare SDOF directly through the modal path by seeding an initial
  // condition: run modal_dynamic with zero load but a nonzero start is not exposed, so
  // instead drive it as a free oscillator via the ModalSystem's SDOF closed form by a
  // brief impulse-equivalent. Simpler: check the transfer via a static-then-release is
  // out of scope; validate the PERIOD through the zero-crossing of a cosine response to
  // a suddenly-applied constant load (step response oscillates about the static value
  // with period 2π/ω).
  numerics::ModalLoad load;
  load.pattern = {1.0};  // unit modal load on the single free DOF
  const Real period = 2.0 * M_PI / omega;
  const Real dt = period / 200.0;
  const Real t_end = 5.0 * period;
  std::vector<numerics::ModalTimePoint> hist =
      numerics::modal_dynamic(sys, load, dt, t_end);

  // Undamped step response: u(t) = (F/k)(1 - cos ω t). Static value F/k; the response
  // oscillates between 0 and 2F/k with period 2π/ω. Modal load p = φᵀ pattern; physical
  // u = φ q with q(t) = (p/ω²)(1 - cos ω t), so u = φ p (1-cos)/ω² and φ²/... — compare
  // against the direct closed form at every sample.
  const Real p = basis.mode_free[0][0] * 1.0;  // φ_1 · pattern
  const Real phi = basis.mode_free[0][0];
  Real max_err = 0.0;
  for (const numerics::ModalTimePoint& f : hist) {
    const Real t = f.time;
    const Real q = (p / (omega * omega)) * (1.0 - std::cos(omega * t));
    const Real u_exact = phi * q;
    max_err = std::max(max_err, std::fabs(f.displacement[0][0] - u_exact));
  }
  // Exact integrator: agreement to machine-ish precision over 5 periods.
  CX_CHECK(max_err < 1e-9);
}

// (4.1) Damped SDOF step response matches the closed-form under-damped solution
// u(t) = (F/k)[1 - e^{-ζωt}(cos ω_d t + ζω/ω_d sin ω_d t)].
void test_damped_step_response() {
  const Real k = 100.0, mass = 1.0;  // ω = 10 rad/s
  numerics::EigenBasis basis = make_sdof(k, mass);
  const Real omega = basis.modes[0].omega;
  const Real zeta = 0.05;
  numerics::Damping damp;
  damp.modal_ratios = {zeta};
  numerics::ModalSystem sys = numerics::project_modal_system(basis, damp);
  CX_NEAR(sys.zeta[0], zeta, 1e-12);

  numerics::ModalLoad load;
  load.pattern = {1.0};
  const Real period = 2.0 * M_PI / omega;
  const Real dt = period / 300.0;
  const Real t_end = 8.0 * period;
  std::vector<numerics::ModalTimePoint> hist =
      numerics::modal_dynamic(sys, load, dt, t_end);

  const Real phi = basis.mode_free[0][0];
  const Real p = phi;  // φ·1
  const Real wd = omega * std::sqrt(1.0 - zeta * zeta);
  Real max_err = 0.0, u_static = 0.0;
  for (const numerics::ModalTimePoint& f : hist) {
    const Real t = f.time;
    const Real q_static = p / (omega * omega);
    const Real q = q_static *
                   (1.0 - std::exp(-zeta * omega * t) *
                              (std::cos(wd * t) +
                               zeta * omega / wd * std::sin(wd * t)));
    const Real u_exact = phi * q;
    u_static = phi * q_static;
    max_err = std::max(max_err, std::fabs(f.displacement[0][0] - u_exact));
  }
  CX_CHECK(max_err < 1e-8 * std::fabs(u_static) + 1e-12);
  // Late-time response settles toward the static deflection within the remaining decay
  // envelope e^{-ζω t_end} (≈ 0.08 here) — the transient has not fully died at 8 periods.
  const Real envelope = std::exp(-zeta * omega * t_end);
  CX_CHECK(std::fabs(hist.back().displacement[0][0] - u_static) <=
           1.1 * envelope * std::fabs(u_static));
}

// (4.2) Steady-state harmonic response: static deflection at low Ω, resonant peak
// F/k · Q with Q = 1/(2ζ), and half-power bandwidth Δω ≈ 2ζω.
void test_steady_state_resonance() {
  const Real k = 100.0, mass = 1.0;  // ω = 10 rad/s
  numerics::EigenBasis basis = make_sdof(k, mass);
  const Real omega = basis.modes[0].omega;
  const Real zeta = 0.02;
  numerics::Damping damp;
  damp.modal_ratios = {zeta};
  numerics::ModalSystem sys = numerics::project_modal_system(basis, damp);

  const std::vector<Real> pattern = {1.0};  // unit force on the free DOF
  const Real phi = basis.mode_free[0][0];

  // Static (Ω -> 0): |u| -> φ²·|pattern|/ω² = F/k (single mode, phi²=1/mass=1 here).
  numerics::HarmonicResponse dc =
      numerics::steady_state_response(sys, pattern, 1e-4);
  const Real u_static = std::abs(dc.amplitude[0]);
  CX_NEAR(u_static, phi * phi / (omega * omega), 1e-6);

  // At resonance Ω = ω: |u| = static · Q, Q = 1/(2ζ).
  numerics::HarmonicResponse res =
      numerics::steady_state_response(sys, pattern, omega);
  const Real u_peak = std::abs(res.amplitude[0]);
  const Real Q = 1.0 / (2.0 * zeta);
  CX_NEAR(u_peak, u_static * Q, 1e-4 * u_static * Q);

  // Half-power bandwidth: the two frequencies where |u| = u_peak/sqrt(2) straddle ω by
  // Δω ≈ 2ζω (so ω_hi - ω_lo ≈ 2ζω). Locate them by a fine scan.
  const Real target = u_peak / std::sqrt(2.0);
  Real w_lo = 0.0, w_hi = 0.0;
  Real prev_w = 0.5 * omega, prev_a = std::abs(
      numerics::steady_state_response(sys, pattern, prev_w).amplitude[0]);
  const int N = 4000;
  for (int i = 1; i <= N; ++i) {
    const Real w = 0.5 * omega + (1.5 * omega - 0.5 * omega) * i / N;
    const Real a =
        std::abs(numerics::steady_state_response(sys, pattern, w).amplitude[0]);
    if (w_lo == 0.0 && prev_a < target && a >= target)
      w_lo = prev_w + (w - prev_w) * (target - prev_a) / (a - prev_a);
    if (w_lo != 0.0 && w_hi == 0.0 && prev_a >= target && a < target)
      w_hi = prev_w + (w - prev_w) * (target - prev_a) / (a - prev_a);
    prev_w = w;
    prev_a = a;
  }
  CX_CHECK(w_lo > 0.0 && w_hi > w_lo);
  const Real bandwidth = w_hi - w_lo;
  CX_NEAR(bandwidth, 2.0 * zeta * omega, 0.05 * 2.0 * zeta * omega);
}

// (1.5) Modal load projection p_k = φ_kᵀ f and (4.3) base-motion effective load
// p_k = -Γ_k, on a 2-DOF spring-mass chain (two modes).
void test_projection_and_base_motion() {
  const Real k = 100.0, mass = 2.0;
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
  m.procedure = Procedure::Frequency;

  const fem::LinearSystem K = fem::assemble_linear_static(m);
  const fem::LinearSystem M = fem::assemble_mass(m, false);
  numerics::EigenBasis basis = numerics::extract_modes(K, M, 2);
  numerics::Damping damp;
  numerics::ModalSystem sys = numerics::project_modal_system(basis, damp);

  // Project an arbitrary free-DOF load and compare to the direct dot product φ_kᵀ f.
  const std::vector<Real> f_free = {3.0, -1.5};
  const std::vector<Real> p = sys.project_load(f_free);
  for (std::size_t kk = 0; kk < 2; ++kk) {
    Real dot = 0.0;
    for (std::size_t i = 0; i < 2; ++i) dot += basis.mode_free[kk][i] * f_free[i];
    CX_NEAR(p[kk], dot, 1e-12);
  }

  // Base motion in x: effective load f_eff = -M r; projected p_k = -Γ_k (participation).
  const std::vector<Real> f_base = numerics::base_motion_load(basis, M, 0);
  const std::vector<Real> p_base = sys.project_load(f_base);
  const numerics::Participation part = numerics::participation(basis, M, 0);
  for (std::size_t kk = 0; kk < 2; ++kk)
    CX_NEAR(p_base[kk], -part.factor[kk], 1e-9);
}

// (4.3) Rayleigh damping C = αM + βK maps to modal ζ_k = (α/ω_k + βω_k)/2.
void test_rayleigh_mapping() {
  numerics::Damping d;
  d.alpha = 0.7;
  d.beta = 1.3e-4;
  const Real w = 250.0;
  CX_NEAR(d.ratio(0, w), 0.5 * (d.alpha / w + d.beta * w), 1e-14);
  // Explicit modal ratios override Rayleigh.
  d.modal_ratios = {0.03};
  CX_NEAR(d.ratio(0, w), 0.03, 1e-14);
  CX_NEAR(d.ratio(1, w), 0.5 * (d.alpha / w + d.beta * w), 1e-14);  // beyond list
}

}  // namespace

int main() {
  test_free_period_undamped();
  test_damped_step_response();
  test_steady_state_resonance();
  test_projection_and_base_motion();
  test_rayleigh_mapping();
  CX_MAIN_RETURN();
}
