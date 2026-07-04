# Design Sensitivity and Optimization

## MODIFIED Requirements

### Requirement: Design variables and responses
The solver SHALL define coordinate design variables via `*DESIGN VARIABLES, TYPE=COORDINATE`, a named design response via `*DESIGN RESPONSE` / `*OBJECTIVE` (compliance / strain energy `fᵀu`, or a single nodal displacement DOF `cᵀu`), and SHALL associate the response with the design variables for sensitivity evaluation, exposing that association through the Python bindings. Orientation design variables, filters, feasible-direction optimization, and robust/random-field design are deferred and, when present, the deck is rejected with a clear message rather than silently mis-handled. (ref: src/designvariabless.f, src/objectives.f)

#### Scenario: Declare a coordinate design response
- GIVEN `*DESIGN VARIABLES, TYPE=COORDINATE` naming a node set and a `*DESIGN RESPONSE` (strain energy or a nodal displacement)
- WHEN the model is read
- THEN each named node's x, y, z coordinate SHALL become a design variable
- AND the response SHALL be associated with those design variables
- AND the design-variable list and response SHALL be inspectable from the Python bindings

#### Scenario: Deferred design feature is rejected, not faked
- GIVEN a deck with `*DESIGN VARIABLES, TYPE=ORIENTATION`, `*FILTER`, `*FEASIBLE DIRECTION`, or `*ROBUST DESIGN`
- WHEN the model is read
- THEN the solver SHALL reject the deck with a message naming the unsupported card
- AND SHALL NOT report a fabricated sensitivity result

### Requirement: Sensitivity analysis
A `*SENSITIVITY` step SHALL compute the gradient of the defined response with respect to each coordinate design variable by the adjoint method — reusing the primal linear-static solution (and, for a general linear response, one adjoint solve with the SAME assembled operator) — via NumPP/SciPP through the ComputeBackend, and the reported gradient SHALL match a central finite-difference of the response (perturb the coordinate by ±h, re-solve, `(O(+h)−O(−h))/2h`) to better than 1e-4 relative on a small linear-elastic model. (ref: src/sensitivitys.f, src/objective_shapeener_tot.f)

#### Scenario: Adjoint gradient matches finite difference
- GIVEN a `*SENSITIVITY` step with a strain-energy or nodal-displacement `*DESIGN RESPONSE` and coordinate design variables on a small cantilever beam
- WHEN the step runs
- THEN the solver SHALL output dObjective/dx for each design variable
- AND each dObjective/dx SHALL agree with the central finite-difference gradient of the response to within 1e-4 relative

#### Scenario: Adjoint reuses the primal operator
- GIVEN a nodal-displacement response with N coordinate design variables
- WHEN the sensitivity step runs
- THEN the primal system SHALL be assembled and solved once
- AND the adjoint response solve SHALL reuse that same assembled operator (not re-factor per design variable)
