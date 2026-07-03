#pragma once
#include <complex>
#include <cstddef>
#include <functional>
#include <vector>

#include "calculixpp/core/types.hpp"
#include "calculixpp/fem/assembly.hpp"
#include "calculixpp/numerics/eigensolution.hpp"

// Modal superposition engine (spec: eigensolution — modal superposition projection;
// dynamic-analysis — *MODAL DYNAMIC / *STEADY STATE DYNAMICS / damping / base motion).
//
// A mass-normalized eigenbasis Φ (φᵀ M φ = 1, φᵀ K φ = λ = ω²) decouples the equations
// of motion. Projecting u = Φ q turns  M ü + C u̇ + K u = f(t)  into N independent SDOF
// modal equations
//     q̈_k + 2 ζ_k ω_k q̇_k + ω_k² q_k = p_k(t),   p_k(t) = φ_kᵀ f(t),
// under proportional (Rayleigh or modal) damping. *MODAL DYNAMIC integrates each SDOF
// in time and recombines u = Σ q_k φ_k; *STEADY STATE DYNAMICS solves the same SDOFs in
// the frequency domain (complex amplitude per mode) over a frequency sweep. This is the
// reusable modal machinery layered on top of the eigensolution basis — no re-factorization
// of the full system. (ref: src/dyna.c, src/steadystate.c)
namespace cxpp::numerics {

// Damping model for the modal SDOFs. Rayleigh damping C = alpha*M + beta*K maps to a
// per-mode ratio ζ_k = (alpha/ω_k + beta*ω_k) / 2 (spec: dynamic-analysis — Rayleigh
// damping). Explicit per-mode modal damping ratios (*MODAL DAMPING) override this when
// supplied (modal_ratios non-empty, indexed by mode). If both are zero/empty the system
// is undamped.
struct Damping {
  Real alpha{0.0};                  // Rayleigh mass coefficient (C += alpha M)
  Real beta{0.0};                   // Rayleigh stiffness coefficient (C += beta K)
  std::vector<Real> modal_ratios;   // explicit ζ_k per mode (overrides Rayleigh if set)

  // Effective modal damping ratio ζ_k for mode k with angular frequency omega_k.
  Real ratio(std::size_t k, Real omega_k) const;
};

// The reduced modal system for a mass-normalized basis: per-mode angular frequency
// ω_k, modal stiffness ω_k², unit modal mass (mass-normalized), and the modal damping
// ratio ζ_k. This is the projection of the operators onto modal coordinates (task 1.5).
// `project_load` turns a physical free-DOF load into the modal load vector p_k = φ_kᵀ f.
struct ModalSystem {
  std::size_t n_modes{0};
  std::vector<Real> omega;   // ω_k (rad/time)
  std::vector<Real> lambda;  // ω_k²  (modal stiffness, unit modal mass)
  std::vector<Real> zeta;    // ζ_k  (modal damping ratio)
  const EigenBasis* basis{nullptr};  // borrowed; the mode_free vectors project loads

  // Project a physical free-DOF load vector (length basis->n_free) onto the modal
  // coordinates: p_k = φ_kᵀ f. Returns one modal load per mode.
  std::vector<Real> project_load(const std::vector<Real>& free_load) const;
};

// Build the reduced modal system from a mass-normalized eigenbasis and a damping model
// (task 1.5 — modal superposition projection of operators). `basis` must outlive the
// returned ModalSystem (it is borrowed for load projection).
ModalSystem project_modal_system(const EigenBasis& basis, const Damping& damping);

// ---------------------------------------------------------------------------
// *MODAL DYNAMIC — transient response by modal superposition (task 4.1).
// ---------------------------------------------------------------------------

// A physical nodal load pattern scaled by a time function. `pattern` is the full
// free-DOF load vector (length n_free); the applied load at time t is
// pattern * amplitude(t). A step (constant) load uses amplitude ≡ 1.
struct ModalLoad {
  std::vector<Real> pattern;  // free-DOF spatial load pattern (constant in time)
};

// One modal-dynamic output frame: the physical nodal displacement at a sampled time.
struct ModalTimePoint {
  Real time{0.0};
  std::vector<Vec3> displacement;  // full nodal displacement (size num_nodes)
};

// Integrate the decoupled modal SDOFs of `sys` under the spatial load `load` scaled by
// the piecewise-linear time function `amplitude` sampled at the step times, from t=0 to
// t_end with fixed step `dt`, recombining the physical displacement at each step. Each
// SDOF is advanced with the EXACT closed-form recurrence for a linearly-varying load over
// the step (the analytical impulse-response / Duhamel solution), so the integrator is
// exact for any load that is piecewise-linear in time — the natural period is reproduced
// with no numerical-damping error, which is what makes the analytical SDOF validation
// exact. Zero initial displacement/velocity (rest start). (spec: dynamic-analysis —
// *MODAL DYNAMIC.)
//
// `amplitude(t)` returns the load scale at time t (default: constant 1 — a step load).
std::vector<ModalTimePoint> modal_dynamic(
    const ModalSystem& sys, const ModalLoad& load, Real dt, Real t_end,
    const std::function<Real(Real)>& amplitude = nullptr);

// ---------------------------------------------------------------------------
// *STEADY STATE DYNAMICS — harmonic response over a frequency sweep (task 4.2).
// ---------------------------------------------------------------------------

// Harmonic response at one excitation frequency: the complex physical displacement
// amplitude per node (magnitude = peak amplitude, arg = phase lead) for a load
// f(t) = pattern * cos(Ω t). Recombined from the per-mode complex modal amplitudes
// q_k = p_k / (ω_k² - Ω² + 2 i ζ_k ω_k Ω).
struct HarmonicResponse {
  Real frequency{0.0};                        // excitation frequency Ω/(2π) (cycles/time)
  Real omega{0.0};                            // excitation Ω (rad/time)
  std::vector<std::complex<Real>> amplitude;  // per free-DOF complex amplitude
};

// Compute the steady-state harmonic response of `sys` to the spatial load `pattern`
// (free-DOF, applied as pattern * cos(Ω t)) at excitation angular frequency `omega_exc`.
// Returns the complex free-DOF amplitude; magnitude is the response amplitude and phase
// is arg. Resonance peaks occur at ω_exc → ω_k, limited by the modal damping ζ_k (the
// peak amplitude of an undamped-driven mode is 1/(2ζ_k) times its static response).
// (spec: dynamic-analysis — steady-state dynamics.)
HarmonicResponse steady_state_response(const ModalSystem& sys,
                                       const std::vector<Real>& pattern,
                                       Real omega_exc);

// Sweep the steady-state response across `num` excitation frequencies logarithmically
// (num >= 2) or a single point (num == 1) between f_lo and f_hi (cycles/time). Returns
// one HarmonicResponse per frequency. (spec: dynamic-analysis — frequency sweep.)
std::vector<HarmonicResponse> steady_state_sweep(const ModalSystem& sys,
                                                 const std::vector<Real>& pattern,
                                                 Real f_lo, Real f_hi, std::size_t num);

// ---------------------------------------------------------------------------
// *BASE MOTION — support excitation (task 4.3).
// ---------------------------------------------------------------------------

// Build the effective free-DOF load pattern for a rigid-base acceleration of unit
// magnitude along translational direction `dir` (0=x,1=y,2=z): f_eff = -M r, where r is
// the rigid-body influence vector for that direction. Projected onto mode k this gives
// the modal participation load p_k = -Γ_k (Γ_k the participation factor), so a base
// acceleration a(t) drives each SDOF by -Γ_k a(t). Used by *BASE MOTION to convert a
// prescribed support acceleration into a modal forcing. (spec: dynamic-analysis — base
// motion.) `M` must be the same mass system used to build `basis`.
std::vector<Real> base_motion_load(const EigenBasis& basis,
                                   const fem::LinearSystem& M, int dir);

}  // namespace cxpp::numerics
