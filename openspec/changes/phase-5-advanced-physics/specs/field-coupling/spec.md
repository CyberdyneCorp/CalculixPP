## ADDED Requirements

### Requirement: Segregated multi-field iteration
The engine SHALL orchestrate a segregated solve of coupled physical fields, iterating each field to its own convergence and exchanging coupling terms between fields until the coupled (outer) residual converges. Each field solve reuses nonlinear-solution-control and linear-algebra-and-solvers. (ref: src/compfluidfem.c, src/electromagnetics.c)

#### Scenario: Outer coupling loop converges
- GIVEN two coupled fields with an initial state
- WHEN the segregated coupling loop runs
- THEN each field SHALL be solved in turn with updated coupling terms until the outer coupled residual falls below the coupling tolerance

#### Scenario: Under-relaxation on a stiff coupling
- GIVEN a strongly coupled field pair prone to oscillation
- WHEN under-relaxation is configured
- THEN the exchanged coupling terms SHALL be relaxed between outer iterations to stabilize convergence

### Requirement: Coupling-term exchange
The engine SHALL transfer coupling quantities between fields at their shared interfaces or overlapping domains — for example Joule heat from the electromagnetic field to the thermal field, temperature to temperature-dependent fluid/material properties, forced-convection film between a network and a solid, and fluid tractions to a structure. (ref: src/radflowload.c, src/electromagnetics.c)

#### Scenario: Joule heat drives the thermal field
- GIVEN an electromagnetic solve producing ohmic (Joule) losses
- WHEN the thermal field is updated
- THEN the Joule heat SHALL be applied as a volumetric heat source to the thermal field (see electromagnetic-analysis, heat-transfer-analysis)

### Requirement: Pressure–velocity coupling for CFD
The engine SHALL provide the pressure-velocity coupling (SIMPLE-class) for the incompressible and compressible Navier-Stokes fields, with each pressure/velocity solve executed through NumPP and the ComputeBackend. (ref: src/compfluidfem.c)

#### Scenario: Incompressible pressure correction
- GIVEN a momentum predictor step for an incompressible flow
- WHEN the pressure-velocity coupling runs
- THEN a pressure-correction solve SHALL enforce continuity and update the velocity field

### Requirement: Coupled convergence control
The engine SHALL judge convergence of the coupled system across all participating fields using an outer coupling tolerance in addition to each field's inner convergence, and SHALL cut back or relax when the coupled residual stalls. (ref: src/compfluidfem.c)

#### Scenario: Coupled residual stall triggers relaxation
- GIVEN a coupled solve whose outer residual stops decreasing
- WHEN the stall is detected
- THEN the engine SHALL increase under-relaxation or reduce the increment rather than declare false convergence

### Requirement: Conjugate interface consistency
The engine SHALL enforce flux and state continuity (e.g. heat flux and temperature) at conjugate fluid–solid and network–solid interfaces, sharing one interface geometry between the coupled fields. (ref: src/radflowload.c)

#### Scenario: Conjugate heat-transfer interface
- GIVEN a fluid and a solid sharing a conjugate interface
- WHEN the coupled step converges
- THEN heat flux and interface temperature SHALL be continuous across the interface within the coupling tolerance

### Requirement: Reuse across coupled procedures
The engine SHALL serve 3-D CFD, electromagnetic-thermal (Joule) coupling, and 1-D network–thermal coupling through one orchestration, so the segregated coupling logic is specified once and reused. (ref: cfd-and-network-analysis, electromagnetic-analysis, heat-transfer-analysis)

#### Scenario: Same engine for EM-thermal and network-thermal
- GIVEN an electromagnetic-thermal deck and, separately, a network-thermal deck
- WHEN each is solved
- THEN both SHALL use the same field-coupling orchestration with their respective coupling terms
