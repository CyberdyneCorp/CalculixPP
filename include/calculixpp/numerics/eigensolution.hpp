#pragma once
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

}  // namespace cxpp::numerics
