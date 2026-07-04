# eigensolution Specification

## Purpose
TBD - created by archiving change phase-4-dynamics-and-eigenproblems. Update Purpose after archive.
## Requirements
### Requirement: Eigenpair extraction and mass normalization
The engine SHALL extract the requested number of generalized symmetric eigenpairs `K x = λ M x` through the linear-algebra-and-solvers NumPP eigensolver, mass-normalize the eigenvectors, and expose them as a modal basis. Extraction uses subspace iteration / Lanczos with a spectral shift, executed through the ComputeBackend. (ref: src/arpack.c, src/arpackcs.c)

#### Scenario: Extract and normalize requested modes
- GIVEN a symmetric `K`, `M` and a requested number of modes
- WHEN extraction runs
- THEN the engine SHALL return that many eigenpairs with mass-normalized eigenvectors (`xᵀ M x = 1`)

#### Scenario: Modes returned in ascending order
- GIVEN a completed extraction
- WHEN the modal basis is exposed
- THEN the modes SHALL be ordered by ascending eigenvalue

### Requirement: Spectral shift and rigid-body handling
The engine SHALL apply a spectral shift so that rigid-body (zero-frequency) and buckling (non-positive) eigenvalues are extracted robustly via shift-invert on a NumPP factorization. (ref: src/arpackbu.c)

#### Scenario: Rigid-body modes extracted
- GIVEN an unconstrained model with rigid-body modes
- WHEN extraction runs with the default shift
- THEN the near-zero-frequency rigid-body modes SHALL be found without failure of the factorization

### Requirement: Modal participation and effective mass
The engine SHALL compute modal participation factors and modal effective mass for a specified excitation direction from the mass-normalized basis. (ref: src/dyna.c)

#### Scenario: Participation for a base excitation
- GIVEN a mass-normalized modal basis and an excitation direction
- WHEN participation is requested
- THEN the engine SHALL return per-mode participation factors and effective masses whose sum approaches the total mass as modes are added

### Requirement: Modal superposition projection
The engine SHALL project the system operators and applied loads onto the modal basis to form the reduced modal equations, so modal-dynamic and steady-state-dynamic procedures integrate in modal coordinates. (ref: src/dyna.c, src/steadystate.c)

#### Scenario: Project loads onto modal coordinates
- GIVEN a modal basis and a physical load vector
- WHEN the system is projected
- THEN the engine SHALL produce the reduced modal load and diagonal (or block) modal operators for superposition

### Requirement: Complex and damped modes
The engine SHALL solve the complex eigenproblem arising from proportional damping by reducing the damping operator onto a preceding real mass-normalized eigenbasis and solving the resulting quadratic modal eigenproblem, returning complex eigenpairs (damped natural frequency, damping ratio, decay rate) and complex mode shapes for the complex-frequency procedure.

This is an ABAQUS-style subspace-projected complex-modes capability for proportional (Rayleigh / modal) damping. It is deliberately NOT the CalculiX `*COMPLEX FREQUENCY, CORIOLIS` gyroscopic problem (which reduces a skew operator and solves `(-ω²I + Λ + i·ω·G_r)q = 0`); that path is a documented follow-on. The reduced quadratic problem `(λ²I + λC_r + Λ)q = 0` is linearized to a real `2·nev × 2·nev` companion pencil and solved with the dense NumPP `eig()` non-symmetric eigensolver on the small reduced matrix — no sparse complex eigensolver is required. (ref: NumPP linalg.hpp `eig`; CalculiX complexfreq.c / coriolissolve.f as the contrasting gyroscopic path)

#### Scenario: Complex eigenpairs from proportional damping
- GIVEN a real mass-normalized eigenbasis (`φᵀ M φ = 1`) and a proportional (Rayleigh or modal) damping model
- WHEN a complex-frequency extraction runs
- THEN the engine SHALL return complex eigenvalues `λ` whose negative real part gives the decay rate and whose imaginary part gives the damped angular frequency, together with complex mode shapes recombined from the real basis (`φ_c = Φ q`)

#### Scenario: Rayleigh damping gives the diagonal closed form
- GIVEN a Rayleigh damping `C = α M + β K`
- WHEN the reduced modal problem is formed
- THEN the reduced modal damping SHALL be diagonal (`C_r,kk = α + β·ω_k²`) and each mode's complex eigenvalue SHALL equal the closed form `λ_k = -ζ_k·ω_k ± i·ω_k·√(1-ζ_k²)` with `ζ_k = (α/ω_k + β·ω_k)/2` within the documented tolerance

#### Scenario: Conjugate-pair de-duplication and ordering
- GIVEN the requested number of complex modes
- WHEN eigenvalues are returned
- THEN they SHALL be de-duplicated to one representative per complex-conjugate pair (`Im(λ) ≥ 0`), overdamped/real roots (`ω_d = 0`) SHALL be reported rather than dropped, and the modes SHALL be ordered by ascending undamped magnitude `|λ|`

### Requirement: Shared basis across procedures
The engine SHALL provide the same extracted real eigenbasis to the frequency, buckling, complex-frequency, modal-dynamic, steady-state-dynamic, and substructure consumers so eigenextraction is performed once per configuration and reused; the complex-frequency consumer SHALL reuse the preceding `*FREQUENCY` extraction to build the reduced complex problem rather than re-extracting or performing a full-scale complex eigensolve. (ref: modal-and-buckling-analysis, dynamic-analysis, substructure-generation)

#### Scenario: Reuse a basis for modal dynamics
- GIVEN a completed `*FREQUENCY` extraction preceding a `*MODAL DYNAMIC` step
- WHEN the modal-dynamic step runs
- THEN it SHALL reuse the existing modal basis rather than re-extracting the eigenpairs

#### Scenario: Reuse the real basis for the complex-frequency reduction
- GIVEN a `*FREQUENCY` step preceding a `*COMPLEX FREQUENCY` step
- WHEN the complex-frequency step runs
- THEN it SHALL reuse the existing real modal basis (`Φ`) to build the reduced complex problem and SHALL NOT re-extract the real eigenpairs or perform a full-scale complex eigensolve

