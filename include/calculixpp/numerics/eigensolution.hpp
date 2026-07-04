#pragma once
#include <complex>
#include <cstddef>
#include <optional>
#include <vector>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/types.hpp"
#include "calculixpp/fem/assembly.hpp"

// Eigensolution engine (spec: eigensolution [NEW]). The reusable modal machinery on
// top of the dense generalized symmetric eigensolve: extract the lowest N eigenpairs
// of K x = λ M x, mass-normalize them (xᵀ M x = 1), sort ascending, and expose
// participation factors / modal effective mass. Consumed by *FREQUENCY (and later by
// buckling, complex-frequency, modal/steady-state dynamics, and substructure). The
// raw eigensolve routes through NumPP (numpp::linalg::eigh / cholesky) — the DENSE
// generalized path is used for correctness+validation; a sparse shift-invert Lanczos
// on SciPP's factorization is the scalable target (see the .cpp note).
namespace cxpp::numerics {

// The proportional-damping model (Rayleigh α/β + explicit modal ζ_k) that the complex
// modal reduction projects onto the real basis. Defined in modal_dynamics.hpp; forward
// declared here to avoid a circular include (that header includes this one).
struct Damping;

// One extracted mode: its eigenvalue λ (= ω² for a *FREQUENCY problem), the derived
// natural angular frequency ω = sqrt(λ) (rad/time) and cyclic frequency f = ω/(2π)
// (cycles/time), and the FULL nodal mode shape (mass-normalized, size num_nodes),
// with constrained DOFs expanded through the constraint transform (SPC value 0 in a
// modal problem, MPC slaves reconstructed from masters).
struct Mode {
  Real eigenvalue{0.0};   // λ  (ω² for free vibration)
  Real omega{0.0};        // sqrt(max(λ,0))  — angular frequency (rad/time)
  Real frequency{0.0};    // ω / (2π)         — cyclic frequency (cycles/time)
  std::vector<Vec3> shape;  // mass-normalized nodal displacement mode shape
};

// Result of a generalized eigenextraction: the requested modes (ascending λ) plus
// the free-DOF numbering used, so consumers (modal dynamics, participation) can
// project physical vectors onto the modal basis. `mode_free[k]` is mode k's
// eigenvector on the free DOFs (size n_free), mass-normalized (φᵀ M φ = 1).
struct EigenBasis {
  std::vector<Mode> modes;               // ascending eigenvalue
  Index n_free{0};                       // free-DOF count of the K/M system
  std::vector<std::vector<Real>> mode_free;  // per-mode free-DOF eigenvector (φ)
  fem::LinearSystem stiffness;           // the assembled K system (dof_eq/transform)
};

// Extract the lowest `num_modes` generalized eigenpairs of the assembled stiffness
// `K` and mass `M` (both on the SAME free-DOF numbering — pass the systems from
// assemble_linear_static / assemble_mass for one model). Dense path: Cholesky
// M = L Lᵀ, reduce to the standard symmetric A = L⁻¹ K L⁻ᵀ, eigh(A), back-transform
// φ = L⁻ᵀ z, mass-normalize, sort ascending. A spectral `shift` σ (default 0) solves
// (K - σ M) instead so rigid-body / near-zero modes are handled robustly; the
// returned eigenvalues are un-shifted (λ = μ + σ). Throws std::runtime_error if M is
// not positive definite (Cholesky fails) or num_modes exceeds the free-DOF count.
EigenBasis extract_modes(const fem::LinearSystem& K, const fem::LinearSystem& M,
                         std::size_t num_modes, Real shift = 0.0);

// Convenience: assemble K and M from `model` and extract `num_modes`. Uses the
// consistent mass matrix. (spec: eigensolution — shared basis; modal-and-buckling —
// *FREQUENCY.)
EigenBasis extract_modes(const Model& model, std::size_t num_modes, Real shift = 0.0);

// Extract the lowest `num_modes` linear-buckling factors of the pencil
// (K + λ K_geo) φ = 0 (spec: modal-and-buckling-analysis — *BUCKLE). This is a
// DEDICATED extractor, NOT extract_modes: the buckling eigenvalue is a load factor λ,
// not ω², so there is no ≥0 clamp and no omega/frequency mapping — the returned Mode's
// `eigenvalue` field holds λ (ascending POSITIVE), `omega`/`frequency` are left 0, and
// `shape` is the buckling mode.
//
// The pencil is solved as (A = −K_geo, B = K) with B = K the SPD UNLOADED stiffness
// (the classical linear-buckling anchor): Cholesky K = L Lᵀ, reduce to the standard
// symmetric Â = L⁻¹(−K_geo)L⁻ᵀ, eigh(Â) → θ, back-transform φ = L⁻ᵀ z, map λ = 1/θ
// (since (K + λ K_geo)φ = 0 ⟺ (−K_geo)φ = (1/λ) K φ, so θ = 1/λ). θ ≤ 0 (rigid-body /
// near-null K_geo directions / load reversal) maps to λ ≤ 0 or λ → +∞ and is rejected;
// only positive λ survive, sorted ascending. `K` and `Kgeo` must share the free-DOF
// numbering (assemble both from one model). The DENSE generalized path is used; the
// scalable sparse path is gated on upstream SciPP target-selection (SciPP#18).
EigenBasis extract_buckling_modes(const fem::LinearSystem& K,
                                  const fem::LinearSystem& Kgeo,
                                  std::size_t num_modes);

// Per-mode modal participation factors and modal effective mass for a rigid-body
// excitation direction (spec: eigensolution — participation and effective mass). For
// the mass-normalized mode φ_k and the unit rigid-body direction vector d (a nodal
// field with `dir` = 1 on one translational DOF of every node), the participation
// factor is Γ_k = φ_kᵀ M r, where r is the influence vector (the free-DOF restriction
// of d), and the modal effective mass is m_eff_k = Γ_k². The effective masses sum to
// the total translational mass in that direction as the basis is completed.
struct Participation {
  int direction{0};                 // excitation DOF direction 0..2 (x,y,z)
  std::vector<Real> factor;         // Γ_k per mode
  std::vector<Real> effective_mass; // m_eff_k = Γ_k² per mode
  Real total_effective_mass{0.0};   // Σ_k m_eff_k (approaches the rigid-body mass)
};

// Compute participation factors / effective mass of `basis` for excitation along
// translational direction `dir` (0=x,1=y,2=z), using the mass system `M`. (Requires
// the same M passed to extract_modes.)
Participation participation(const EigenBasis& basis, const fem::LinearSystem& M,
                            int dir);

// ---------------------------------------------------------------------------
// Damped complex modes (*COMPLEX FREQUENCY, proportional damping — option B).
// ---------------------------------------------------------------------------

// One damped complex mode of (λ²M + λC + K)x = 0. For proportional (Rayleigh / modal)
// damping the reduced problem is diagonal, so each mode has the exact SDOF closed form
// λ = -ζω_n ± i ω_n √(1-ζ²). We keep the representative with Im(λ) ≥ 0 of each
// conjugate pair. `omega_d` = |Im(λ)| (damped angular frequency), `omega_n` = |λ|
// (undamped angular frequency), `zeta` = -Re(λ)/|λ| (>0 stable/decaying, <0 growing),
// `frequency` = ω_d/(2π), `decay_rate` = Re(λ). The physical complex mode shape is
// φ_c = Φ q (upper reduced block back to physical), split into real/imag full nodal
// fields.
struct ComplexMode {
  std::complex<Real> eigenvalue{0.0, 0.0};  // λ
  Real omega_d{0.0};                        // |Im(λ)| — damped angular frequency
  Real omega_n{0.0};                        // |λ|     — undamped angular frequency
  Real zeta{0.0};                           // -Re(λ)/|λ| — damping ratio
  Real frequency{0.0};                      // ω_d / (2π) — damped cyclic frequency
  Real decay_rate{0.0};                     // Re(λ)   — real part (decay if <0)
  std::vector<Vec3> shape_real;             // Re(φ_c) full nodal field
  std::vector<Vec3> shape_imag;             // Im(φ_c) full nodal field
};

// Result of the damped complex-mode reduction: the complex modes (ascending |λ|) and
// the free-DOF count of the underlying system.
struct ComplexEigenBasis {
  std::vector<ComplexMode> modes;
  Index n_free{0};
};

// Reduce the proportional damping operator onto the real mass-normalized eigenbasis
// `real_basis` (Φᵀ M Φ = I, Φᵀ K Φ = Λ = diag(ω_k²)) and extract the lowest
// `num_modes` damped complex modes. Forms the reduced quadratic (λ²I + λC_r + Λ)q = 0
// with C_r = Φᵀ C Φ (diagonal α + β·ω_k² for Rayleigh, overridden by 2·ζ_k·ω_k where
// modal ratios are set), linearizes to the real 2·nev companion A = [[0,I],[-Λ,-C_r]],
// and solves it with numpp::linalg::eig. This is the proportional-damping option-(B)
// path: it is EXACT for proportional damping (diagonal C_r) and is explicitly NOT the
// CalculiX CORIOLIS gyroscopic problem (a skew G_r with an i·ω coupling — a different
// eigenproblem). The reduced-operator assembly carries an (empty here) skew/imaginary
// block so a future gyroscopic G_r plugs into the same linearization without redesign.
ComplexEigenBasis extract_complex_modes(const EigenBasis& real_basis,
                                        const Damping& damping, std::size_t num_modes);

// Guarded small-problem dense oracle for cross-validation only: solve the full physical
// 2n state-space companion [[0,I],[-M⁻¹K,-M⁻¹C]] with numpp::linalg::eig, where
// C = αM + βK is the Rayleigh damping of `damping` (modal ratios are ignored — this
// oracle validates the Rayleigh reduction). Throws if n_free exceeds kDenseComplexMaxDof
// or if the modal-ratio path is requested. Returns the same post-processed complex modes
// (ascending |λ|, Im(λ) ≥ 0) so the modal-reduced result can be compared entry-for-entry.
ComplexEigenBasis extract_complex_modes_dense(const fem::LinearSystem& K,
                                              const fem::LinearSystem& M,
                                              const Damping& damping,
                                              std::size_t num_modes);

}  // namespace cxpp::numerics
