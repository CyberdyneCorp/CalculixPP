# Electromagnetic Analysis

## Purpose

CalculixPP solves electromagnetic field problems — magnetostatics, electromagnetics
(eddy currents), and electric conduction — ported to pure C++20 and selected by
`*ELECTROMAGNETICS`, optionally coupled to heat transfer through Joule heating. The
potential-field linear/nonlinear solves run through **NumPP** (replacing the reference
sparse solvers) with assembly/solve dispatched via the **ComputeBackend** (CPU always
available; GPU optional). (ref: src/electromagnetics.c, src/electromagneticss.f,
src/resultsem.f, src/magneticpermeabilitys.f, src/electricalconductivitys.f)

**Porting Phase:** 5

## Requirements

### Requirement: Electromagnetic procedure
An `*ELECTROMAGNETICS` step SHALL solve the electromagnetic field for the magnetic
and/or electric potentials, supporting magnetostatic, harmonic, and transient
sub-types. The field system SHALL be assembled and solved through NumPP on the
ComputeBackend, and the step SHALL be reachable from the Python bindings.

#### Scenario: Magnetostatic solve
- GIVEN an `*ELECTROMAGNETICS` step of magnetostatic type with coil currents
- WHEN the step runs
- THEN CalculixPP SHALL output the magnetic field / flux density distribution

### Requirement: Electromagnetic material properties
Electromagnetic analysis SHALL consume `*MAGNETIC PERMEABILITY` and
`*ELECTRICAL CONDUCTIVITY` material properties, and related constants where
applicable, matching the reference card semantics.

#### Scenario: Permeability drives the field
- GIVEN materials with differing `*MAGNETIC PERMEABILITY`
- WHEN the magnetostatic field is solved
- THEN the flux density distribution SHALL reflect the per-material permeabilities

### Requirement: Coupling to thermal field
Electromagnetic analysis SHALL support coupling to heat transfer so that resistive
(Joule) heating from the electric field acts as a heat source in the thermal solution
(see heat-transfer-analysis).

#### Scenario: Joule heating
- GIVEN an electromagnetic step coupled to heat transfer with an electric current
- WHEN the step runs
- THEN resistive heating SHALL contribute to the temperature field

### Requirement: Air/surrounding domain
The formulation SHALL allow an enveloping air domain so that the magnetic field can
be computed in the region surrounding conductors and magnetic bodies. The air-domain
mesh SHALL be prepared via CyberCadKernel (see mesh-processing).

#### Scenario: Field in surrounding air
- GIVEN a conductor enclosed by an air domain
- WHEN the magnetostatic field is solved
- THEN the magnetic field SHALL be computed in the air region surrounding the conductor

### Requirement: Reference fidelity
Electromagnetic results SHALL match the reference CalculiX output for the
corresponding `test/` deck within the documented numerical tolerance, on the CPU
ComputeBackend with no GPU present and on any optional GPU backend.

#### Scenario: Regression against reference deck
- GIVEN a reference electromagnetic `*.inp` deck with known CalculiX results
- WHEN CalculixPP solves it on the CPU backend
- THEN the field / flux density distribution SHALL agree with the reference within tolerance
