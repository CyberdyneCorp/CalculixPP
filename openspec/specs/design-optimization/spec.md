# Design Sensitivity and Optimization

## Purpose

Computes design sensitivities of objectives and constraints with respect to design variables
(coordinates or orientations), drives gradient-based design improvement, and supports robust
design under uncertainty. In CalculixPP the sensitivity/adjoint linear solves and filtering
operators are pure C++20 kernels on **NumPP** dense/sparse containers, assembled and solved
through the **ComputeBackend** (CPU default and always available; GPU optional and never
required). All optimization operations are reachable from the Python bindings.
(ref: src/CalculiX.c:1555 sensi_coor / sensi_orien, src/CalculiX.c:1607 robustdesign,
src/CalculiX.c:1642 feasibledirection, src/sensitivitys.f, src/objectives.f,
src/designvariabless.f)

**Porting Phase:** 5

## Requirements

### Requirement: Design variables and responses
The solver SHALL define design variables via `*DESIGN VARIABLES` (coordinate-based shape
variables or orientation variables), design responses via `*DESIGN RESPONSE`, objectives via
`*OBJECTIVE`, and constraints via `*CONSTRAINT` / `*GEOMETRIC CONSTRAINT`, preserving the
reference keyword semantics. (ref: src/designvariabless.f, src/objectives.f)

#### Scenario: Declare an objective
- GIVEN `*DESIGN VARIABLES` and an `*OBJECTIVE` (e.g. mass or compliance)
- WHEN the model is read
- THEN the objective SHALL be associated with the design variables for sensitivity evaluation
- AND the association SHALL be inspectable from the Python bindings

### Requirement: Sensitivity analysis
A `*SENSITIVITY` step SHALL compute the gradients of the defined objective(s) and
constraint(s) with respect to the design variables, using coordinate sensitivities
(`sensi_coor`) or orientation sensitivities (`sensi_orien`) as appropriate, with the adjoint
linear solves performed via NumPP through the ComputeBackend. (ref: src/sensitivitys.f)

#### Scenario: Compute objective gradient
- GIVEN a `*SENSITIVITY` step with a defined objective and design variables
- WHEN the step runs
- THEN the solver SHALL output the gradient of the objective with respect to each design variable

### Requirement: Filtering and geometric tolerances
Sensitivity/optimization SHALL support `*FILTER` (sensitivity filtering over a radius) and
`*GEOMETRIC TOLERANCE` to regularize and constrain the design field, implemented as C++20
forward/backward filter operators on NumPP containers. (ref: src/filters.f,
src/filterbackwardmain.c, src/filterforwardmain.c)

#### Scenario: Sensitivity filtering
- GIVEN a `*FILTER` with a radius applied to a sensitivity field
- WHEN sensitivities are computed
- THEN they SHALL be smoothed over the filter radius before use

### Requirement: Feasible direction optimization
A `*FEASIBLE DIRECTION` step SHALL compute an improved design direction that decreases the
objective while respecting active constraints, using the sensitivity information.
(ref: src/feasibledirection.c)

#### Scenario: Improved design direction
- GIVEN computed objective and constraint sensitivities
- WHEN a `*FEASIBLE DIRECTION` step runs
- THEN the solver SHALL output a design-update direction that reduces the objective while respecting active constraints

### Requirement: Robust design / random field
A `*ROBUST DESIGN` step SHALL evaluate the sensitivity of the response to random input
variation defined via `*RANDOM FIELD` and `*CORRELATION LENGTH`, supporting uncertainty-aware
design, with the random-field decomposition computed via NumPP. (ref: src/robustdesign.c,
src/randomfieldmain.c)

#### Scenario: Response under random variation
- GIVEN a `*ROBUST DESIGN` step with a `*RANDOM FIELD` and `*CORRELATION LENGTH`
- WHEN the step runs
- THEN the solver SHALL report the sensitivity of the response to the modeled input uncertainty

### Requirement: Reference-result fidelity
Sensitivity gradients, filtered fields, and feasible-direction updates for a given reference `test/` deck SHALL match the reference CalculiX output within the documented numerical tolerance, on the CPU backend with no GPU present. (ref: src/sensitivitys.f)

#### Scenario: Match reference sensitivity deck
- GIVEN a reference sensitivity/optimization `*.inp` deck run on the CPU ComputeBackend with no GPU toolkit installed
- WHEN the analysis completes
- THEN the computed gradients and design updates SHALL agree with the reference CalculiX results within tolerance
