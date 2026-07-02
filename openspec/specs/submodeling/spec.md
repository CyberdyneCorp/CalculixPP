# Submodeling

## Purpose

CalculixPP supports submodeling — driving a local, finer model from the results of a
prior global analysis at the submodel's cut boundary — ported to pure C++20. A
`*SUBMODEL` card declares the global boundary (nodes and/or surfaces) whose driven
values are read from the global solution and interpolated onto the submodel boundary;
`*BOUNDARY, SUBMODEL` and `*CLOAD/*DSLOAD, SUBMODEL` then apply those driven values as
the local model's boundary conditions (see loads-and-boundary-conditions). Interpolation
and the local solve run through the `ComputeBackend` (CPU/NumPP by default, GPU optional
and never required). Submodel results SHALL match reference CalculiX within tolerance and
are reachable from the Python bindings. (ref: src/submodels.f, src/calinput.f)

**Porting Phase:** 5 — Advanced physics

## Requirements

### Requirement: Submodel boundary declaration
A `*SUBMODEL` card SHALL declare the submodel's cut boundary as a set of nodes and/or surfaces and identify the global result source (the global `.frd`/results file and the global element set searched for host elements). The declaration SHALL be reachable from the Python bindings. (ref: src/submodels.f, src/calinput.f)

#### Scenario: Declare node-based submodel boundary
- GIVEN a `*SUBMODEL, TYPE=NODE` card naming a boundary node set and a global result file
- WHEN the deck is parsed
- THEN the solver SHALL register those nodes as the driven boundary sourced from that global solution

#### Scenario: Declare surface-based submodel boundary
- GIVEN a `*SUBMODEL, TYPE=SURFACE` card naming a boundary surface and a global result file
- WHEN the deck is parsed
- THEN the solver SHALL register that surface as the driven boundary sourced from that global solution

### Requirement: Interpolation of global results onto the submodel boundary
Driven values at each submodel boundary node SHALL be obtained by locating the host global element containing the node and interpolating the global solution field with that element's shape functions. (ref: src/submodels.f)

#### Scenario: Interpolate driven displacements
- GIVEN a submodel boundary node lying inside a global element
- WHEN driven values are computed
- THEN the boundary node's driven displacement SHALL equal the global displacement field evaluated at the node's location via the host element shape functions

#### Scenario: Node outside the global element set
- GIVEN a submodel boundary node that falls outside every element in the referenced global element set
- WHEN driven values are computed
- THEN the program SHALL report an error that no host element was found for that node

### Requirement: Driven boundary conditions and loads
`*BOUNDARY, SUBMODEL` SHALL apply the interpolated global displacements as prescribed DOF values on the submodel boundary, and `*CLOAD, SUBMODEL` / `*DSLOAD, SUBMODEL` SHALL apply interpolated global forces and pressures, consistent with loads-and-boundary-conditions. (ref: src/submodels.f, src/calinput.f)

#### Scenario: Prescribed submodel displacements
- GIVEN a `*BOUNDARY, SUBMODEL` card referencing a declared submodel node set and a global step
- WHEN the submodel step runs
- THEN the listed DOFs SHALL be prescribed to the interpolated global displacement values

#### Scenario: Driven distributed load
- GIVEN a `*DSLOAD, SUBMODEL` card on a declared submodel surface
- WHEN the submodel step runs
- THEN the surface SHALL be loaded with the interpolated global pressure field

### Requirement: Submodel solve
With its boundary driven from the global results, the submodel SHALL run as an ordinary analysis step (static or other supported procedure) on the finer local mesh, with assembly and solve routed through the `ComputeBackend`. (ref: src/submodels.f)

#### Scenario: No GPU present
- GIVEN a build with no CUDA/OpenCL/Metal backend available and a driven submodel step
- WHEN the step runs
- THEN the submodel SHALL be assembled and solved on the CPU/NumPP backend and produce correct local results

### Requirement: Reference-result fidelity
For every submodel `test/` deck, CalculixPP submodel boundary values and local results SHALL match reference CalculiX output within a documented relative tolerance, validated through the pytest bindings.

#### Scenario: Submodel results match reference
- GIVEN a reference CalculiX submodel deck and its expected local results
- WHEN CalculixPP runs the same deck driven by the same global solution
- THEN the computed submodel results SHALL match the reference within the documented relative tolerance
