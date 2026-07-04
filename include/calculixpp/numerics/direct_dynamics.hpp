#pragma once
#include <functional>
#include <optional>
#include <vector>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/types.hpp"
#include "calculixpp/fem/assembly.hpp"
#include "calculixpp/numerics/linear_static.hpp"
#include "calculixpp/numerics/modal_dynamics.hpp"

// Direct-integration dynamics (spec: dynamic-analysis — Direct-integration dynamics;
// tasks 3.1 / 3.2 / 3.3). Integrates the full equations of motion
//     M a + C v + K u = f(t)
// in physical coordinates by the implicit HHT-α (Hilber-Hughes-Taylor) generalization
// of Newmark-β, time-stepped by fixed increments. Unlike modal superposition this keeps
// the full (sparse) operators and re-solves each step, so it carries nonlinear material
// / contact tangents and inertia together (task 3.2). Assembly and the per-step linear
// solve route through fem::assemble_* and numerics::solve_reduced (ComputeBackend /
// SciPP); no GPU is required. (ref: src/nonlingeo.c, src/dynamics.f)
//
// HHT-α enforces equilibrium at t_{n+1+α}:
//     M a_{n+1} + (1+α)(C v_{n+1} + K u_{n+1}) - α(C v_n + K u_n)
//         = (1+α) f_{n+1} - α f_n,
// with the Newmark update  u_{n+1} = u_n + Δt v_n + Δt²[(½-β) a_n + β a_{n+1}],
//                          v_{n+1} = v_n + Δt[(1-γ) a_n + γ a_{n+1}].
// For controllable numerical damping α ∈ [-1/3, 0]; β = (1-α)²/4, γ = ½ - α give
// unconditional stability and second-order accuracy. α = 0 is the trapezoidal
// (average-acceleration) Newmark rule — no algorithmic damping, so an undamped free
// vibration conserves energy over the run (task 3.3, the energy check).
namespace cxpp::numerics {

// HHT-α integration parameters. `alpha` is the numerical-damping knob (α ≤ 0; α = 0 is
// energy-conserving trapezoidal Newmark, more negative α damps high-frequency response
// more strongly). β and γ are derived from α unless overridden. (task 3.3 — numerical
// damping control.)
struct HhtParams {
  Real alpha{0.0};    // HHT α ∈ [-1/3, 0]
  Real beta{0.25};    // Newmark β (default (1-α)²/4, filled by make_hht)
  Real gamma{0.5};    // Newmark γ (default ½-α, filled by make_hht)

  // Newmark β = (1-α)²/4, γ = ½-α for the given α. Clamps α to [-1/3, 0].
  static HhtParams from_alpha(Real alpha);
};

// One output frame of a direct-dynamics run: physical nodal displacement / velocity /
// acceleration at a sampled time, plus the energy balance at that time (task 3.3).
struct DirectTimePoint {
  Real time{0.0};
  std::vector<Vec3> displacement;  // full nodal displacement (size num_nodes)
  std::vector<Vec3> velocity;      // full nodal velocity
  std::vector<Vec3> acceleration;  // full nodal acceleration
  Real kinetic_energy{0.0};        // ½ vᵀ M v (free DOFs)
  Real strain_energy{0.0};         // ½ uᵀ K u (linear) / internal strain energy
  Real total_energy{0.0};          // kinetic + strain
};

// Diagnostics from a direct-dynamics run (steps taken, Newton iterations for the
// nonlinear path, and the peak relative drift of total energy over an undamped run —
// the energy-conservation measure used by task 3.3).
struct DirectReport {
  int steps{0};
  int iterations{0};          // total Newton iterations (nonlinear path); 0 if linear
  bool nonlinear{false};
  Real energy_drift{0.0};     // max |E(t) - E(0)| / E(0) over the run (0 if E(0)==0)
};

// Options controlling a direct-dynamics run.
struct DirectOptions {
  HhtParams hht{};             // integration scheme (default α=0 trapezoidal Newmark)
  Damping damping{};           // Rayleigh C = αM + βK (modal_ratios unused here)
  bool lumped_mass{false};     // row-sum lumped mass instead of consistent
  bool nonlinear{false};       // include material/contact tangent each Newton iteration
  std::optional<SolverKind> forced{};  // override the linear solver

  // Optional non-rest initial conditions on the FREE DOFs (size n_free; empty = 0).
  // A nonzero initial displacement/velocity with zero load gives a free vibration whose
  // total energy (½vᵀMv + ½uᵀKu) is conserved for the undamped α=0 scheme — the
  // energy-conservation validation of task 3.3.
  std::vector<Real> u0_free;   // initial free-DOF displacement
  std::vector<Real> v0_free;   // initial free-DOF velocity
};

// Integrate a `*DYNAMIC` step of `model` from t=0 to `t_end` with fixed step `dt`,
// starting from rest (u=v=0, a from M a₀ = f(0)). The spatial load pattern is the
// model's external load at full magnitude, scaled in time by `amplitude(t)` (default:
// constant — a suddenly-applied step load at t=0⁺). Returns one DirectTimePoint per
// step (including t=0). With opts.nonlinear the per-step effective tangent adds the
// material-point tangent (and contact) to the inertial/damping terms (task 3.2);
// otherwise the linear K is used and factored implicitly each step through
// solve_reduced. (spec: dynamic-analysis — implicit time step.)
std::vector<DirectTimePoint> direct_dynamic(
    const Model& model, Real dt, Real t_end, const DirectOptions& opts = {},
    const std::function<Real(Real)>& amplitude = nullptr,
    DirectReport* report = nullptr);

}  // namespace cxpp::numerics
