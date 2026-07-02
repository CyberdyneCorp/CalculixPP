# Dynamic Analysis

## Purpose

CalculixPP computes the time- or frequency-domain response of a structure to dynamic
loading, either by direct time integration or by modal superposition, ported to pure
C++20. Assembly and any linear solves run through the `ComputeBackend` abstraction
(CPU via NumPP by default; GPU optional and never required), and modal procedures
consume eigenmodes produced by the sibling **modal-and-buckling-analysis** capability.
Results SHALL match reference CalculiX for the corresponding `test/` deck within a
documented numerical tolerance. (ref: src/CalculiX.c, src/nonlingeo.c, src/dyna.c, src/steadystate.c)

**Porting Phase:** 4 — Dynamics & eigenproblems

## Requirements

### Requirement: Direct-integration dynamics
A `*DYNAMIC` step SHALL integrate the equations of motion `M a + C v + K u = f` in
time. Implicit integration (alpha-method / HHT) SHALL be the default; explicit
central-difference integration SHALL be selected by `*DYNAMIC, EXPLICIT`. Time-step
assembly and linear solves SHALL be performed through the `ComputeBackend` (CPU/NumPP
by default), and the step SHALL be reachable from the Python bindings. (ref: src/nonlingeo.c, src/dynamics.f)

#### Scenario: Implicit time step
- GIVEN an implicit `*DYNAMIC` step
- WHEN an increment is solved
- THEN inertia and damping contributions SHALL be included
- AND Newton-Raphson iterations SHALL converge the increment before advancing time

#### Scenario: Explicit stability
- GIVEN a `*DYNAMIC, EXPLICIT` step
- WHEN integrating
- THEN the stable time increment SHALL be governed by the element wave-speed / critical-time-step estimate

#### Scenario: No GPU present
- GIVEN a build with no CUDA/OpenCL/Metal backend available
- WHEN a `*DYNAMIC` step runs
- THEN the step SHALL execute on the CPU/NumPP backend and produce correct results

### Requirement: Modal dynamic analysis
A `*MODAL DYNAMIC` step (nmethod 4) SHALL compose the transient response by
superposition of the eigenmodes from a preceding `*FREQUENCY` step (see
modal-and-buckling-analysis), applying modal damping, without re-factorizing the full
system. The modal reduction and integration SHALL be C++20 code using NumPP for any
linear algebra. (ref: src/dyna.c, src/modaldynamics.f)

#### Scenario: Modal transient response
- GIVEN a `*FREQUENCY` step followed by a `*MODAL DYNAMIC` step
- WHEN the dynamic step runs
- THEN the response SHALL be computed from the modal coordinates, not by re-factorizing the full system

#### Scenario: Missing eigenmode source
- GIVEN a `*MODAL DYNAMIC` step with no preceding `*FREQUENCY` step in the same job
- WHEN the job runs
- THEN the program SHALL report an error that the required eigenmodes are unavailable

### Requirement: Steady-state dynamics
A `*STEADY STATE DYNAMICS` step (nmethod 5) SHALL compute the harmonic (frequency-
domain) response across a requested frequency range by modal superposition, with
damping, producing amplitude and phase results. It SHALL require a preceding
`*FREQUENCY` step (see modal-and-buckling-analysis). (ref: src/steadystate.c, src/steadystatedynamicss.f)

#### Scenario: Frequency sweep
- GIVEN a `*STEADY STATE DYNAMICS` step over a frequency range
- WHEN the step runs
- THEN the solver SHALL output the harmonic response amplitude and phase at the requested frequencies

### Requirement: Damping
Dynamic procedures SHALL support Rayleigh damping (`*DAMPING, ALPHA=, BETA=`) and
`*MODAL DAMPING` for modal procedures. (ref: src/dampings.f, src/modaldampings.f)

#### Scenario: Rayleigh damping
- GIVEN a dynamic step with `*DAMPING, ALPHA=, BETA=`
- WHEN the equations of motion are integrated
- THEN the damping matrix SHALL be formed as `alpha*M + beta*K`

### Requirement: Base motion
Modal procedures SHALL support `*BASE MOTION` to apply prescribed support excitation
(e.g. for seismic / shaker analysis). (ref: src/basemotions.f)

#### Scenario: Support excitation
- GIVEN a modal procedure with `*BASE MOTION` prescribing support acceleration
- WHEN the step runs
- THEN the response SHALL be computed relative to the prescribed base excitation

### Requirement: Reference-result fidelity
For every dynamic `test/` deck (direct, modal, steady-state), CalculixPP results SHALL
match reference CalculiX output within a documented numerical tolerance (displacement /
amplitude / phase time histories within relative tolerance), validated through the
pytest bindings.

#### Scenario: Deck matches reference
- GIVEN a reference CalculiX dynamic test deck and its expected results
- WHEN CalculixPP runs the same deck
- THEN the computed response SHALL match the reference within the documented relative tolerance

### Requirement: Green-function step (Phase 4)
A `*GREEN` step SHALL compute the structural response to unit excitations (Green functions) at prescribed degrees of freedom, producing the response basis used for subsequent steady-state / harmonic dynamics. The unit-excitation solves SHALL run through NumPP on the ComputeBackend (CPU by default, GPU optional and never required) and be reachable from the Python bindings. (ref: src/steadystate.c, src/greens.f)

#### Scenario: Response to a unit excitation
- GIVEN a `*GREEN` step prescribing unit excitations at selected degrees of freedom
- WHEN the step runs
- THEN a Green function SHALL be computed for each unit excitation through the ComputeBackend
- AND the resulting response basis SHALL be available to a following steady-state dynamics step

#### Scenario: No GPU present
- GIVEN a build with no CUDA/OpenCL/Metal backend available
- WHEN a `*GREEN` step runs
- THEN the unit-excitation solves SHALL execute on the CPU/NumPP backend and produce correct results
