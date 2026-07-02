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
A `*BUCKLE` step (nmethod 3) SHALL compute buckling load factors and buckling mode
shapes by solving `(K + lambda K_geo) phi = 0` with a NumPP eigensolver, where the
geometric stiffness derives from a reference load in the step. (ref: src/arpackbu.c, src/bucklings.f)

#### Scenario: Buckling factor
- GIVEN a `*BUCKLE` step with a reference load
- WHEN the step runs
- THEN the solver SHALL output buckling load factors; the critical load is the factor times the reference load

### Requirement: Complex frequency analysis
A `*COMPLEX FREQUENCY` step (nmethod 6/7) SHALL compute complex eigenmodes built from
the real eigenmodes of a preceding `*FREQUENCY` step, accounting for damping, Coriolis,
or friction effects, using a NumPP complex eigensolver. (ref: src/complexfreq.c, src/complexfrequencys.f)

#### Scenario: Damped complex modes
- GIVEN a `*FREQUENCY` step followed by `*COMPLEX FREQUENCY` with damping
- WHEN the step runs
- THEN the solver SHALL output complex eigenvalues whose imaginary and real parts give damped frequency and decay rate

### Requirement: Preceding frequency step requirement
Modal-superposition procedures SHALL require a preceding `*FREQUENCY` step in the same
job that produced the eigenmodes — this applies to modal dynamic and steady-state
dynamics (see dynamic-analysis) and to complex frequency. (ref: src/CalculiX.c, src/dyna.c, src/steadystate.c)

#### Scenario: Missing eigenmode source
- GIVEN a `*COMPLEX FREQUENCY` step with no preceding `*FREQUENCY` step in the job
- WHEN the job runs
- THEN the program SHALL report an error that the required eigenmodes are unavailable

### Requirement: Reference-result fidelity
For every eigenvalue `test/` deck (frequency, buckling, complex frequency), CalculixPP eigenfrequencies, buckling factors, and complex eigenvalues SHALL match reference CalculiX output within a documented relative tolerance, validated through the pytest bindings.

#### Scenario: Eigenfrequencies match reference
- GIVEN a reference CalculiX `*FREQUENCY` test deck and its expected eigenfrequencies
- WHEN CalculixPP runs the same deck
- THEN each computed eigenfrequency SHALL match the reference within the documented relative tolerance
