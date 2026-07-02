# Loads and Boundary Conditions

## Purpose

Defines the prescribed displacements/temperatures, concentrated and distributed loads,
thermal loads, and their time variation (amplitudes), applied within an analysis step.
CalculixPP parses the same Abaqus-style cards as reference CalculiX and assembles the
right-hand side and prescribed-DOF constraints for the linear `*STATIC` solve `K u = f`.
Results SHALL match reference CalculiX for the corresponding test deck within documented
tolerance. (ref: src/boundarys.f, src/cloads.f, src/dloads.f, src/amplitudes.f,
src/temperatures.f)

**Porting Phase:** 1 — Foundation & vertical slice (*BOUNDARY, *CLOAD, *DLOAD pressure,
*TEMPERATURE, and constant/step *AMPLITUDE); body loads, thermal fluxes/films/radiation,
and general amplitude curves land in Phase 2 (noted per requirement).

## Requirements

### Requirement: Boundary conditions (prescribed DOF) (Phase 1)
The model SHALL support `*BOUNDARY` to prescribe degrees of freedom (displacements,
rotations, temperatures, or other DOFs) at nodes or node sets, with a value (default
zero) over a DOF range. Boundary conditions given in the model-data section SHALL be
fixed (homogeneous) for all steps unless overridden inside a step. Prescribed DOFs SHALL
be enforced in the linear `*STATIC` solve. (ref: src/boundarys.f)

#### Scenario: Prescribed nonzero displacement
- GIVEN `*BOUNDARY` inside a step with `Nset, 1, 1, 0.5`
- WHEN the step is solved
- THEN DOF 1 of the nodes in Nset SHALL be enforced to 0.5
- AND the resulting displacement field SHALL match reference CalculiX within tolerance

### Requirement: Concentrated loads (Phase 1)
The model SHALL support `*CLOAD` (concentrated force/moment on a DOF) at nodes or node
sets, and `*CFLUX` (concentrated heat flux, Phase 3). Concentrated loads SHALL be added
to the right-hand-side vector `f`. (ref: src/cloads.f, src/cfluxes.f)

#### Scenario: Point force on a node
- GIVEN `*CLOAD` with `node, 2, 100.`
- WHEN the step is solved
- THEN a force of 100 SHALL be applied in DOF 2 of the node
- AND the response SHALL match reference CalculiX within tolerance

### Requirement: Distributed loads (Phase 1: pressure)
The model SHALL support `*DLOAD` face pressure `P` on element faces (Phase 1), and body
loads (gravity, centrifugal), `*DSLOAD`, `*DFLUX`, `*FILM`, and `*RADIATE` (Phase 2/3).
Pressure SHALL be integrated over the loaded element faces into consistent nodal forces
via a C++20 kernel on the CPU ComputeBackend. (ref: src/dloads.f, src/dfluxes.f,
src/films.f, src/radiates.f)

#### Scenario: Pressure on an element face
- GIVEN `*DLOAD` with `elset, P3, 1.0e5`
- WHEN the step is solved
- THEN a pressure of 1.0e5 SHALL be applied to face 3 of the elements in elset
- AND the equivalent nodal forces SHALL match reference CalculiX within tolerance

### Requirement: Thermal loading
The model SHALL support `*TEMPERATURE` to prescribe nodal temperatures (a load in a
mechanical step, a boundary value in a heat-transfer step) and `*INITIAL CONDITIONS,
TYPE=TEMPERATURE` for the reference/initial field. Thermal strain SHALL be computed
against the material `*EXPANSION` (see material-models). (ref: src/temperatures.f)

#### Scenario: Thermal expansion from a temperature load
- GIVEN a mechanical step with `*TEMPERATURE` raising nodal temperatures and a material `*EXPANSION`
- WHEN the step is solved
- THEN thermal strains SHALL be generated relative to the initial temperature field

### Requirement: Amplitudes and time variation (Phase 1: constant/step)
The model SHALL support `*AMPLITUDE` defining a time-dependent multiplier curve,
referenceable by `AMPLITUDE=` on load and boundary cards. At minimum, constant/step
amplitude and linear ramp over a static step SHALL be supported in Phase 1; general
tabular/frequency curves land in Phase 2. By default loads ramp linearly over a static
step and are applied as steps in dynamics. (ref: src/amplitudes.f)

#### Scenario: Load following an amplitude
- GIVEN a `*CLOAD` referencing an `*AMPLITUDE` curve
- WHEN the load is evaluated at a given step time
- THEN its magnitude SHALL equal the base value times the amplitude value at that time

### Requirement: Load accumulation across steps
Within a step, loads SHALL by default replace those of the previous step unless the step
requests otherwise (`OP=MOD` keeps prior loads, `OP=NEW` removes them), matching
reference CalculiX semantics. (ref: src/calinput.f, src/cloads.f)

#### Scenario: Remove previous loads
- GIVEN a step whose load card specifies `OP=NEW`
- WHEN the step is solved
- THEN loads defined in previous steps SHALL be removed before the new ones are applied

### Requirement: Property changes between steps or on restart (Phase 2)
The model SHALL support `*CHANGE MATERIAL` (card `*CHANGEMATERIAL`) to modify material data, `*CHANGE PLASTIC` (card `*CHANGEPLASTIC`) to modify plastic hardening data, and `*CHANGE SOLID SECTION` (card `*CHANGESOLIDSECTION`) to reassign a solid-section definition, applied between steps or on restart, so subsequent steps use the updated properties while preserving reference CalculiX semantics. (ref: src/changematerials.f, src/changeplastics.f, src/changesolidsections.f)

#### Scenario: Change material data between steps
- GIVEN a `*CHANGE MATERIAL` card in a later step altering an existing material's data
- WHEN that step is solved
- THEN the updated material data SHALL be used from that step onward
- AND earlier steps SHALL remain unaffected

#### Scenario: Reassign a solid section on restart
- GIVEN a restarted job with `*CHANGE SOLID SECTION` reassigning an element set's solid section
- WHEN the restarted step runs
- THEN the elements SHALL use the newly assigned solid-section definition
