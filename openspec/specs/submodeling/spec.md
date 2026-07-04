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
A `*SUBMODEL, TYPE=NODE` card SHALL declare the submodel's cut boundary as a set of driven nodes and identify the global result source (the global displacement field and the global element set searched for host elements), and the declaration SHALL be reachable from the Python bindings.

The `TYPE=NODE` boundary names a node set; the referenced global element set bounds the search for the
host element of each driven node. The global result source is an in-memory `GlobalSolution` (global
mesh + nodal displacement field) supplied through the Python API in this slice; reading the field from a
global `.frd`/results file on disk is a documented follow-on. `TYPE=SURFACE` boundaries are deferred.
(ref: src/submodels.f, src/calinput.f, src/interpolsubmodel.f)

#### Scenario: Declare node-based submodel boundary
- GIVEN a `*SUBMODEL, TYPE=NODE` card naming a boundary node set and a global element set
- WHEN the deck is parsed
- THEN the solver SHALL register those nodes as the driven boundary sourced from that global solution

#### Scenario: Reachable from Python
- GIVEN a submodel deck and an in-memory global solution
- WHEN `solve_submodel` is called through the Python bindings
- THEN it SHALL return the submodel displacement/stress/reaction result dict

### Requirement: Interpolation of global results onto the submodel boundary
Driven values at each submodel boundary node SHALL be obtained by locating the host global element containing the node and interpolating the global displacement field with that element's shape functions.

The host element is the global element (within the referenced global element set) whose inverse
isoparametric map places the node at natural coordinates inside the reference domain within tolerance;
the node's natural coordinates are found by a Newton iteration solving `x(ξ) = X`. The interpolated
displacement is `Σ_a N_a(ξ) U_a` over the host element's nodes, reusing the existing element
shape-function machinery. (ref: src/submodels.f, src/interpolsubmodel.f, src/attach_3d.f)

#### Scenario: Interpolate driven displacements
- GIVEN a submodel boundary node lying inside a global element
- WHEN driven values are computed
- THEN the boundary node's driven displacement SHALL equal the global displacement field evaluated at the node's location via the host element shape functions

#### Scenario: Node outside the global element set
- GIVEN a submodel boundary node that falls outside every element in the referenced global element set
- WHEN driven values are computed
- THEN the program SHALL report an error that no host element was found for that node

### Requirement: Driven boundary conditions and loads
`*BOUNDARY, SUBMODEL` SHALL apply the interpolated global displacements as prescribed DOF values on the submodel boundary, consistent with loads-and-boundary-conditions.

Each `*BOUNDARY, SUBMODEL` DOF is prescribed to the interpolated global displacement component at the
node's location; the value is filled from the global solution at solve time rather than read from the
card. Force- and pressure-driven `*CLOAD, SUBMODEL` / `*DSLOAD, SUBMODEL` (which interpolate the global
stress field) and temperature-driven submodels are deferred. (ref: src/submodels.f, src/calinput.f)

#### Scenario: Prescribed submodel displacements
- GIVEN a `*BOUNDARY, SUBMODEL` card referencing a declared submodel node set
- WHEN the submodel step runs
- THEN the listed DOFs SHALL be prescribed to the interpolated global displacement values

### Requirement: Submodel solve
With its boundary driven from the global results, the submodel SHALL run as an ordinary static analysis step on the finer local mesh, with assembly and solve routed through the `ComputeBackend`. (ref: src/submodels.f)

#### Scenario: No GPU present
- GIVEN a build with no CUDA/OpenCL/Metal backend available and a driven submodel step
- WHEN the step runs
- THEN the submodel SHALL be assembled and solved on the CPU/NumPP backend and produce correct local results

### Requirement: Reference-result fidelity
A submodel driven by a global solution's displacements on its cut boundary SHALL reproduce the global displacement field inside the submodel region within a documented relative tolerance, validated through the C++ unit test and the Python bindings.

#### Scenario: Submodel matches the global field inside the cut region
- GIVEN a beam solved globally and a sub-region cut from it whose cut-boundary nodes are driven by the interpolated global displacements
- WHEN CalculixPP solves the submodel
- THEN the submodel displacement field SHALL match the global displacement field inside the region to relative L2 error below 1e-3

