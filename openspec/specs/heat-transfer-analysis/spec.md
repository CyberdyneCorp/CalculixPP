# Heat Transfer and Coupled Thermomechanical Analysis

## Purpose

CalculixPP solves thermal conduction/convection/radiation fields and their coupling
to the mechanical response, ported to pure C++20. Selected by `*HEAT TRANSFER`
(thermal only) or `*COUPLED TEMPERATURE-DISPLACEMENT` (coupled). All linear/nonlinear
solves and the cavity-radiation matrix run through **NumPP** (replacing SPOOLES/PARDISO)
with assembly/solve dispatched via the **ComputeBackend** abstraction (CPU always
available; GPU optional). (ref: src/heattransfers.f, src/couptempdisps.f,
src/radmatrix.f, src/radflowload.c, src/nonlingeo.c)

**Porting Phase:** 3

## Requirements

### Requirement: Heat transfer analysis
A `*HEAT TRANSFER` step SHALL solve the thermal field for nodal temperatures using
material `*CONDUCTIVITY` and (for transient) `*SPECIFIC HEAT` and `*DENSITY`, subject
to thermal boundary conditions, `*CFLUX`/`*DFLUX`, `*FILM`, and `*RADIATE` loads.
Steady-state and transient modes SHALL both be supported, with the linear/nonlinear
system solved through NumPP via the ComputeBackend. The step SHALL be reachable from
the Python bindings.

#### Scenario: Steady-state conduction
- GIVEN a `*HEAT TRANSFER, STEADY STATE` step with fixed temperatures and a heat flux
- WHEN the step runs
- THEN CalculixPP SHALL solve the steady conduction field and output nodal temperatures and heat flux

#### Scenario: Transient conduction
- GIVEN a transient `*HEAT TRANSFER` step with `*SPECIFIC HEAT` and `*DENSITY` and time parameters
- WHEN the step is time-integrated
- THEN the nodal temperature history SHALL be produced over the increments

### Requirement: Coupled thermomechanical analysis
A `*COUPLED TEMPERATURE-DISPLACEMENT` step SHALL solve the temperature and
displacement fields together, so thermal expansion and temperature-dependent
properties affect the mechanical solution and mechanical dissipation/contact can
affect the thermal solution. The coupled Newton iteration SHALL assemble and solve
through NumPP on the ComputeBackend, preserving the reference coupling semantics.

#### Scenario: Coupled solve
- GIVEN a `*COUPLED TEMPERATURE-DISPLACEMENT` step with a heat flux and a mechanical load
- WHEN an increment is solved
- THEN both the temperature and displacement fields SHALL be updated together with their mutual coupling

### Requirement: Radiation with view factors
Radiation loads SHALL compute or read view factors between radiating surfaces;
`*RADIATE` with cavity radiation SHALL assemble a radiation matrix solved through
NumPP, and view factors MAY be cached/read via `*VIEWFACTOR`. Cavity/surface geometry
preparation SHALL reuse the meshing pipeline (see mesh-processing, CyberCadKernel).

#### Scenario: Cavity radiation
- GIVEN `*RADIATE` surfaces forming a cavity
- WHEN the thermal step runs
- THEN view factors between the surfaces SHALL be computed and used to exchange radiative heat

### Requirement: Convective film coupling to networks
Forced-convection `*FILM` conditions SHALL be able to reference network (fluid)
elements so that the film temperature is taken from the coupled gas/liquid network
solution (see cfd-and-network-analysis). The solid thermal field and the network
SHALL be solved consistently within the coupled step.

#### Scenario: Conjugate heat with a network
- GIVEN a coupled step with `*FILM` referencing a fluid network
- WHEN the step runs
- THEN the solid thermal field and the network flow/temperature SHALL be solved consistently

### Requirement: Reference fidelity
Thermal and coupled results SHALL match the reference CalculiX output for the
corresponding `test/` deck within the documented numerical tolerance, on the CPU
ComputeBackend with no GPU present and on any optional GPU backend.

#### Scenario: Regression against reference deck
- GIVEN a reference thermal `*.inp` deck with known CalculiX results
- WHEN CalculixPP solves it on the CPU backend
- THEN the nodal temperatures and fluxes SHALL agree with the reference within tolerance
