## ADDED Requirements

### Requirement: Newton-Raphson iteration driver
The solver SHALL provide a reusable Newton-Raphson driver that, each increment, assembles the tangent and residual, solves the linearized system through the linear-algebra-and-solvers / ComputeBackend spine, updates the solution, and repeats until the convergence controls are satisfied. The driver is procedure-agnostic and shared by nonlinear statics (Phase 2), thermal/coupled analysis (Phase 3), and nonlinear dynamics (Phase 4). (ref: src/nonlingeo.c)

#### Scenario: Iterate to convergence
- GIVEN a nonlinear increment with an initial guess
- WHEN the driver runs
- THEN it SHALL repeat tangent-assemble → residual → solve → update until the residual satisfies the convergence controls, then accept the increment

#### Scenario: Linear problem is a single iteration
- GIVEN a linear model driven through the nonlinear driver
- WHEN one increment is solved
- THEN the result SHALL equal the direct linear solve (see linear-algebra-and-solvers) within numerical tolerance in a single iteration

### Requirement: Convergence controls
The driver SHALL judge convergence using force-residual and displacement-correction criteria whose tolerances and iteration limits are configurable via `*CONTROLS`, defaulting to documented values when the card is absent. (ref: src/controlss.f)

#### Scenario: Within tolerance accepts
- GIVEN an iteration whose force residual and displacement correction are below the configured tolerances
- WHEN convergence is tested
- THEN the driver SHALL accept the increment and advance

#### Scenario: Exceed iteration limit triggers cutback
- GIVEN an increment that does not converge within the configured maximum iterations
- WHEN the limit is reached
- THEN the driver SHALL discard the increment and request a cutback

### Requirement: Automatic incrementation and cutback
The driver SHALL support automatic time incrementation — growing the increment after easy convergence and cutting it back after non-convergence — as well as `DIRECT` fixed increments and `*TIME POINTS`, and SHALL abort the step if the increment falls below the configured minimum. (ref: src/nonlingeo.c, src/controlss.f)

#### Scenario: Cutback on divergence
- GIVEN automatic incrementation
- WHEN an increment fails to converge
- THEN the increment SHALL be restarted at a smaller size, and the step SHALL abort if the size falls below the configured minimum

#### Scenario: Fixed increment with DIRECT
- GIVEN a step specifying `DIRECT` incrementation
- WHEN the step runs
- THEN the driver SHALL use the fixed user increment without automatic resizing

### Requirement: Amplitude-driven time stepping
The driver SHALL sample `*AMPLITUDE` definitions at each increment so that time-varying loads and boundary conditions are applied at the correct fraction of step time. (ref: src/amplitudes.f)

#### Scenario: Tabular amplitude sampled mid-step
- GIVEN a load referencing a tabular `*AMPLITUDE`
- WHEN an increment lands between two amplitude table points
- THEN the applied load factor SHALL be the interpolated amplitude value at that time

### Requirement: Consistent tangent contract
Element and material kernels SHALL supply the algorithmically consistent tangent to the driver so Newton-Raphson retains quadratic convergence, and a kernel MAY fall back to a numerical tangent that is flagged as such. (ref: src/e_c3d.f, src/umat.f)

#### Scenario: Consistent tangent yields quadratic convergence
- GIVEN a nonlinear material supplying its consistent tangent
- WHEN an increment is iterated
- THEN the residual SHALL decrease at the expected quadratic rate near the solution

### Requirement: Optional line search
The driver SHALL optionally apply a line search that scales the Newton update to improve robustness on strongly nonlinear increments, enabled per analysis and off by default.

#### Scenario: Line search recovers a difficult increment
- GIVEN line search is enabled and a full Newton step would increase the residual
- WHEN the update is computed
- THEN the driver SHALL scale the step to reduce the residual rather than diverge
