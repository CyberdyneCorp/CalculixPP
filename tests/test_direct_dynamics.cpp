// Direct-integration dynamics (Phase 4, tasks 3.1 / 3.2 / 3.3). HHT-α implicit
// integration of M a + C v + K u = f(t). Validated ANALYTICALLY against closed-form
// SDOF / spring-mass results:
//   * (3.1) undamped free vibration reproduces the natural period 2π/ω;
//   * (3.1) a suddenly-applied step load gives u(t) = (F/k)(1 - cos ω t);
//   * (3.3) the undamped α=0 (trapezoidal Newmark) scheme CONSERVES total energy
//     ½vᵀMv + ½uᵀKu over a free-vibration run (the energy check);
//   * (3.3) HHT numerical damping (α<0) DISSIPATES energy over the same free run;
//   * (3.1) Rayleigh damping decays the free-vibration amplitude;
//   * (3.2) the nonlinear path reproduces the linear path on a linear-elastic model.
#include <cmath>
#include <vector>

#include "calculixpp/fem/assembly.hpp"
#include "calculixpp/numerics/direct_dynamics.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

// Single grounded spring-mass on node-1 x (ω = sqrt(k/m)); other DOFs constrained.
// Optional *CLOAD F on node-1 x. Mirrors the modal-dynamics SDOF fixture.
Model make_sdof(Real k, Real mass, Real load = 0.0) {
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
  if (load != 0.0) m.cloads.push_back(Cload{1, 1, load});
  m.procedure = Procedure::Static;  // *DYNAMIC reuses the static assembly
  return m;
}

// Peak |u| of the free x-DOF over a history.
Real peak_ux(const std::vector<numerics::DirectTimePoint>& h) {
  Real p = 0.0;
  for (const auto& tp : h) p = std::max(p, std::fabs(tp.displacement[0][0]));
  return p;
}

// (3.1) Suddenly-applied step load: the exact SDOF response is u(t)=(F/k)(1-cos ωt),
// oscillating between 0 and 2F/k with period 2π/ω. Implicit HHT (α=0) reproduces it to
// the second-order time-discretization error (tightens as dt shrinks).
void test_step_response_undamped() {
  const Real k = 400.0, mass = 1.0, F = 5.0;  // ω = 20 rad/s
  Model m = make_sdof(k, mass, F);
  const Real omega = std::sqrt(k / mass);
  const Real period = 2.0 * M_PI / omega;
  const Real dt = period / 400.0;      // fine step -> tight agreement
  const Real t_end = 4.0 * period;
  numerics::DirectOptions opts;        // α=0, undamped
  const std::vector<numerics::DirectTimePoint> h =
      numerics::direct_dynamic(m, dt, t_end, opts);

  Real max_err = 0.0;
  for (const auto& tp : h) {
    const Real u_exact = (F / k) * (1.0 - std::cos(omega * tp.time));
    max_err = std::max(max_err, std::fabs(tp.displacement[0][0] - u_exact));
  }
  // Second-order accurate: at 400 steps/period the phase/amplitude error is tiny.
  CX_CHECK(max_err < 2e-4);
  // Oscillates about the static value F/k with peak ~2F/k.
  CX_NEAR(peak_ux(h), 2.0 * F / k, 5e-4);
}

// (3.1) Free-vibration period: no load, initial displacement u0, undamped. The response
// is u(t)=u0 cos ωt; find the first zero-up-crossing period from the samples and compare
// to 2π/ω.
void test_free_vibration_period() {
  const Real k = 100.0, mass = 4.0;  // ω = 5 rad/s
  Model m = make_sdof(k, mass);
  const Real omega = std::sqrt(k / mass);
  const Real period = 2.0 * M_PI / omega;
  const Real dt = period / 500.0;
  const Real t_end = 3.0 * period;

  numerics::DirectOptions opts;
  opts.u0_free = {0.01};  // initial displacement on the single free DOF
  const std::vector<numerics::DirectTimePoint> h =
      numerics::direct_dynamic(m, dt, t_end, opts);

  // Analytic u(t) = u0 cos ωt; compare directly.
  Real max_err = 0.0;
  for (const auto& tp : h) {
    const Real u_exact = 0.01 * std::cos(omega * tp.time);
    max_err = std::max(max_err, std::fabs(tp.displacement[0][0] - u_exact));
  }
  CX_CHECK(max_err < 5e-6);
}

// (3.3) ENERGY CONSERVATION: undamped α=0 free vibration keeps total energy constant.
void test_energy_conservation_undamped() {
  const Real k = 250.0, mass = 2.0;
  Model m = make_sdof(k, mass);
  const Real omega = std::sqrt(k / mass);
  const Real period = 2.0 * M_PI / omega;
  const Real dt = period / 100.0;
  const Real t_end = 20.0 * period;  // long run: drift would accumulate if present

  numerics::DirectOptions opts;         // α=0 trapezoidal Newmark
  opts.u0_free = {0.02};
  numerics::DirectReport rep;
  const std::vector<numerics::DirectTimePoint> h =
      numerics::direct_dynamic(m, dt, t_end, opts, nullptr, &rep);

  // Total energy at t=0 is ½ k u0² (all strain energy).
  const Real e0 = h.front().total_energy;
  CX_NEAR(e0, 0.5 * k * 0.02 * 0.02, 1e-12);
  // Trapezoidal Newmark conserves total energy to round-off over 20 periods.
  CX_CHECK(rep.energy_drift < 1e-10);
}

// (3.3) HHT numerical damping (α<0) DISSIPATES energy over the same free run — the
// controllable-damping knob. Energy at the end is strictly below the start.
void test_hht_numerical_damping() {
  const Real k = 250.0, mass = 2.0;
  Model m = make_sdof(k, mass);
  const Real omega = std::sqrt(k / mass);
  const Real period = 2.0 * M_PI / omega;
  const Real dt = period / 40.0;     // coarser step so high-freq damping bites
  const Real t_end = 40.0 * period;

  numerics::DirectOptions opts;
  opts.hht = numerics::HhtParams::from_alpha(-0.3);  // strong numerical damping
  opts.u0_free = {0.02};
  numerics::DirectReport rep;
  const std::vector<numerics::DirectTimePoint> h =
      numerics::direct_dynamic(m, dt, t_end, opts, nullptr, &rep);

  CX_NEAR(opts.hht.beta, (1.3 * 1.3) / 4.0, 1e-12);  // β = (1-α)²/4
  CX_NEAR(opts.hht.gamma, 0.5 + 0.3, 1e-12);         // γ = ½-α
  // Energy decayed measurably (α<0 is dissipative) — drift here is a real energy loss.
  CX_CHECK(h.back().total_energy < 0.99 * h.front().total_energy);
}

// (3.1) Rayleigh damping decays the free-vibration amplitude monotonically (per cycle).
void test_rayleigh_damping_decay() {
  const Real k = 100.0, mass = 1.0;  // ω = 10 rad/s
  Model m = make_sdof(k, mass);
  const Real omega = std::sqrt(k / mass);
  const Real period = 2.0 * M_PI / omega;
  const Real dt = period / 200.0;
  const Real t_end = 6.0 * period;

  numerics::DirectOptions opts;
  opts.damping.beta = 0.02;  // stiffness-proportional: ζ = β ω / 2 = 0.1
  opts.u0_free = {0.05};
  const std::vector<numerics::DirectTimePoint> h =
      numerics::direct_dynamic(m, dt, t_end, opts);

  // Peak amplitude over the first quarter vs. the last quarter: strictly decaying.
  Real first = 0.0, last = 0.0;
  const std::size_t q = h.size() / 4;
  for (std::size_t i = 0; i < q; ++i)
    first = std::max(first, std::fabs(h[i].displacement[0][0]));
  for (std::size_t i = 3 * q; i < h.size(); ++i)
    last = std::max(last, std::fabs(h[i].displacement[0][0]));
  CX_CHECK(last < 0.7 * first);  // ζ=0.1 over ~4-6 cycles loses well over 30%
}

// (3.2) Nonlinear path reproduces the linear path on a linear-elastic model (the driver
// contract: nonlinear-reproduces-linear). Same step-load SDOF, opts.nonlinear=true.
void test_nonlinear_reproduces_linear() {
  const Real k = 400.0, mass = 1.0, F = 5.0;
  Model m = make_sdof(k, mass, F);
  const Real omega = std::sqrt(k / mass);
  const Real period = 2.0 * M_PI / omega;
  const Real dt = period / 200.0;
  const Real t_end = 3.0 * period;

  numerics::DirectOptions lin;
  const std::vector<numerics::DirectTimePoint> hl =
      numerics::direct_dynamic(m, dt, t_end, lin);
  numerics::DirectOptions nl;
  nl.nonlinear = true;
  numerics::DirectReport rep;
  const std::vector<numerics::DirectTimePoint> hn =
      numerics::direct_dynamic(m, dt, t_end, nl, nullptr, &rep);

  CX_CHECK(rep.nonlinear);
  CX_CHECK(hl.size() == hn.size());
  Real max_diff = 0.0;
  for (std::size_t i = 0; i < hl.size(); ++i)
    max_diff = std::max(max_diff,
                        std::fabs(hl[i].displacement[0][0] - hn[i].displacement[0][0]));
  CX_CHECK(max_diff < 1e-9);  // Newton on a linear system matches the direct solve
}

}  // namespace

int main() {
  test_step_response_undamped();
  test_free_vibration_period();
  test_energy_conservation_undamped();
  test_hht_numerical_damping();
  test_rayleigh_damping_decay();
  test_nonlinear_reproduces_linear();
  CX_MAIN_RETURN();
}
