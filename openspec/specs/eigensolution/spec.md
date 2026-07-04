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
The engine SHALL support the complex eigenproblem arising from damping or friction-induced instability, producing complex eigenpairs used by the complex-frequency procedure. (ref: src/complexfreq.c)

#### Scenario: Complex eigenpairs with damping
- GIVEN a system with damping or friction contributions
- WHEN a complex-frequency extraction runs
- THEN the engine SHALL return complex eigenvalues and eigenvectors

### Requirement: Shared basis across procedures
The engine SHALL provide the same extracted basis to the frequency, buckling, complex-frequency, modal-dynamic, steady-state-dynamic, and substructure (fixed-interface modes) consumers, so eigenextraction is performed once per configuration and reused. (ref: modal-and-buckling-analysis, dynamic-analysis, substructure-generation)

#### Scenario: Reuse a basis for modal dynamics
- GIVEN a completed `*FREQUENCY` extraction preceding a `*MODAL DYNAMIC` step
- WHEN the modal-dynamic step runs
- THEN it SHALL reuse the existing modal basis rather than re-extracting the eigenpairs

