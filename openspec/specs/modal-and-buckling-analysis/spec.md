# Eigenvalue Analysis (Frequency, Buckling, Complex Frequency)

## Purpose

CalculixPP computes eigenpairs of the model — natural frequencies and mode shapes,
linear buckling factors, and complex eigenmodes (damping, Coriolis, or friction-induced
instability) — ported to pure C++20. Eigenproblems are solved with **NumPP
eigensolvers** (not ARPACK); assembly of the stiffness, mass, and geometric-stiffness
operators and any linear solves run through the `ComputeBackend` (CPU/NumPP by default,
GPU optional and never required). Computed eigenfrequencies and modes SHALL match
reference CalculiX for the corresponding `test/` deck within tolerance, and feed the
sibling **dynamic-analysis** modal procedures. (ref: src/CalculiX.c, src/arpack.c, src/arpackbu.c, src/complexfreq.c)

**Porting Phase:** 4 — Dynamics & eigenproblems
## Requirements
### Requirement: Eigenfrequency analysis
A `*FREQUENCY` step (nmethod 2) SHALL compute the requested number of natural
frequencies and mode shapes by solving the generalized eigenproblem
`(K - omega^2 M) phi = 0` using a NumPP eigensolver. The eigenmodes SHALL be available
to subsequent modal-dynamic, steady-state, and complex-frequency steps (see
dynamic-analysis), and the step SHALL be reachable from the Python bindings. (ref: src/arpack.c, src/frequencys.f)

#### Scenario: Request N modes
- GIVEN a `*FREQUENCY` step requesting N eigenvalues
- WHEN the step runs
- THEN the solver SHALL output up to N eigenfrequencies and their mode shapes

#### Scenario: No GPU present
- GIVEN a build with no CUDA/OpenCL/Metal backend available
- WHEN a `*FREQUENCY` step runs
- THEN the eigenproblem SHALL be assembled and solved on the CPU/NumPP backend and produce correct results

### Requirement: Cyclic-symmetric eigenfrequencies
When a `*CYCLIC SYMMETRY MODEL` is present, the `*FREQUENCY` step SHALL solve the
complex cyclic-symmetric eigenproblem per selected nodal diameter using a NumPP complex
eigensolver. (ref: src/arpackcs.c)

#### Scenario: Per nodal-diameter modes
- GIVEN a cyclic-symmetric `*FREQUENCY` step with `*SELECT CYCLIC SYMMETRY MODES`
- WHEN the step runs
- THEN eigenfrequencies SHALL be reported for each selected nodal diameter

### Requirement: Linear buckling analysis
A `*BUCKLE` step SHALL compute buckling load factors and buckling mode shapes by solving `(K + λ K_geo) φ = 0`, where `K_geo` derives from a reference load's prestress stress state, returning the buckling factors in ascending positive order with their mode shapes, reachable from the CLI and Python bindings.

The eigenproblem is solved as the pencil `(A = −K_geo, B = K)` with `B = K` SPD (the unloaded
stiffness), mapping `λ = −1/θ` and filtering to positive factors; the critical load is the smallest
positive factor times the reference load. The dense generalized path (NumPP Cholesky of `K` + `eigh`)
is the primary validatable path; the scalable sparse path is gated on upstream SciPP target-selection
support. (ref: src/arpackbu.c, src/bucklings.f, src/e_c3d.f)

#### Scenario: Buckling factors and mode shapes
- GIVEN a `*BUCKLE` step with a reference load
- WHEN the step runs
- THEN the solver SHALL output buckling load factors in ascending positive order with their mode shapes
- AND the critical load SHALL be the smallest positive factor times the reference load

#### Scenario: Matches the beamb reference deck
- GIVEN the stock `beamb` buckling deck (C3D20R Euler column, `*BUCKLE 10, 1.e-2`)
- WHEN CalculixPP runs it
- THEN the reported buckling factors SHALL match `beamb.dat.ref` (λ_1 = 48.15, λ_2 = 106.3) within a documented relative tolerance

#### Scenario: Matches the Euler critical load
- GIVEN a slender column meshed with solid elements under axial compression
- WHEN a `*BUCKLE` step runs
- THEN the lowest factor times the reference load SHALL match the Euler critical load `π²EI/(kL)²` within tolerance for pinned-pinned (`k=1`) and clamped-free (`k=2`)

### Requirement: Complex frequency analysis
A `*COMPLEX FREQUENCY` step SHALL compute complex eigenmodes built from the real eigenmodes of a preceding `*FREQUENCY` step, accounting for proportional (Rayleigh / modal) damping, and report each mode's damped natural frequency, damping ratio, decay rate, and stability.

The shipped physics is proportional-damping subspace-projected complex modes (ABAQUS-style): the real damping operator `C = αM + βK` (and/or `*MODAL DAMPING` ratios) is reduced onto the real modal basis to the quadratic problem `(λ²I + λC_r + Λ)q = 0`, linearized to a real `2·nev` companion pencil, and solved with the dense NumPP `eig()`. This is NOT the CalculiX `*COMPLEX FREQUENCY, CORIOLIS` gyroscopic operator (skew reduced matrix, `i·ω·G_r` coupling, rotor whirl / Campbell split, driven by a rotation body load); the true `CORIOLIS` operator and the `FLUTTER` complex applied-force path are documented follow-ons. (ref: NumPP linalg.hpp `eig`; CalculiX complexfreq.c / coriolissolve.f / complexfrequencys.f as the contrasting gyroscopic path)

#### Scenario: Damped complex modes from a proportional-damping step
- GIVEN a `*FREQUENCY` step followed by `*COMPLEX FREQUENCY` with a proportional `*DAMPING` (Rayleigh) or `*MODAL DAMPING` card
- WHEN the step runs
- THEN the solver SHALL output complex eigenvalues whose imaginary part gives the damped frequency and whose real part gives the decay rate, plus a per-mode damping ratio and a stability flag (`stable` when `ζ > 0`)

#### Scenario: Missing eigenmode source
- GIVEN a `*COMPLEX FREQUENCY` step with no preceding `*FREQUENCY` step in the job
- WHEN the job runs
- THEN the program SHALL report an error that the required eigenmodes are unavailable

### Requirement: Preceding frequency step requirement
Modal-superposition procedures SHALL require a preceding `*FREQUENCY` step in the same
job that produced the eigenmodes — this applies to modal dynamic and steady-state
dynamics (see dynamic-analysis) and to complex frequency. (ref: src/CalculiX.c, src/dyna.c, src/steadystate.c)

#### Scenario: Missing eigenmode source
- GIVEN a `*COMPLEX FREQUENCY` step with no preceding `*FREQUENCY` step in the job
- WHEN the job runs
- THEN the program SHALL report an error that the required eigenmodes are unavailable

### Requirement: Reference-result fidelity
For every real eigenvalue `test/` deck (frequency, buckling), CalculixPP eigenfrequencies and buckling factors SHALL match reference CalculiX output within a documented relative tolerance; for the proportional-damping `*COMPLEX FREQUENCY` procedure, complex eigenvalues (damped frequency and damping ratio) SHALL instead match the analytical closed form within a documented relative tolerance, because CalculiX `*COMPLEX FREQUENCY` solves the different (gyroscopic) eigenproblem and its output would not match. All fidelity checks are validated through the pytest bindings.

#### Scenario: Eigenfrequencies match reference
- GIVEN a reference CalculiX `*FREQUENCY` test deck and its expected eigenfrequencies
- WHEN CalculixPP runs the same deck
- THEN each computed eigenfrequency SHALL match the reference within the documented relative tolerance

#### Scenario: Complex modes match the analytical closed form
- GIVEN a Rayleigh-damped SDOF or small cantilever deck with known `ζ`
- WHEN CalculixPP runs the `*FREQUENCY` + proportional `*COMPLEX FREQUENCY` job
- THEN each computed damped frequency and damping ratio SHALL match the analytical value `λ_k = -ζ_k·ω_k ± i·ω_k·√(1-ζ_k²)`, `ζ_k = (α/ω_k + β·ω_k)/2`, within the documented tolerance, and no CalculiX `*COMPLEX FREQUENCY` reference-deck comparison SHALL be asserted for this procedure

### Requirement: Complex-frequency keyword restrictions
The `*COMPLEX FREQUENCY` parser SHALL reject the `CORIOLIS` and `FLUTTER` keywords with an explicit "not yet implemented" error rather than silently accepting a deck it would mis-solve, because the shipped path implements only proportional-damping complex modes and does not implement the CalculiX gyroscopic Coriolis operator (which requires a rotation body load: rotor speed and axis) or the `FLUTTER` complex applied-force path.

#### Scenario: CORIOLIS keyword rejected
- GIVEN a `*COMPLEX FREQUENCY, CORIOLIS` card
- WHEN the step is parsed
- THEN the parser SHALL reject it with a clear "not yet implemented" message naming the Coriolis gyroscopic operator + rotor-speed/axis input as the missing follow-on, rather than mis-solving with proportional damping

#### Scenario: FLUTTER keyword rejected
- GIVEN a `*COMPLEX FREQUENCY, FLUTTER` card
- WHEN the step is parsed
- THEN the parser SHALL reject it with a clear "not yet implemented" message identifying the complex applied-force path as the missing follow-on

