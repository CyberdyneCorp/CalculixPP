# Static Analysis

## Purpose

Solves the equilibrium of a structure under quasi-static loading in CalculixPP,
either as a single linear solve or an incremental-iterative nonlinear solution,
selected by `*STATIC` inside a step. This is the proven-first vertical slice: the
linear path exercises the full C++20 pipeline (deck parse → assemble via
`ComputeBackend` → solve `K u = f` with NumPP → stress/strain recovery → `.frd`/`.dat`).
Physics and keyword semantics mirror the reference; only the implementation stack
changes. (ref: src/CalculiX.c linstatic/nonlingeo dispatch, ref: src/linstatic.c,
ref: src/nonlingeo.c)

**Porting Phase:** Linear static = Phase 1 (foundation & vertical slice); nonlinear
static, automatic incrementation/cutback, and perturbation = Phase 2.

## Requirements

### Requirement: Linear static analysis
CalculixPP SHALL solve a `*STATIC` step with no nonlinear effects (small displacement,
linear-elastic material, no contact, `NLGEOM` off) by assembling the global linear
stiffness matrix `K` and load vector `f` once and solving `K u = f` in a single
increment using the NumPP sparse direct solver. Assembly SHALL be routed through the
`ComputeBackend` abstraction with the CPU backend as the always-available default,
and the solve SHALL succeed with no GPU toolkit present. (ref: src/linstatic.c)

#### Scenario: Single linear solve
- GIVEN a `*STATIC` step on a linear-elastic model without contact and without `*STEP, NLGEOM`
- WHEN the step runs
- THEN CalculixPP SHALL assemble `K` and `f` once
- AND solve `K u = f` in one increment with the NumPP sparse direct solver
- AND recover element stresses and strains from the displacement solution
- AND write displacements and stresses to the results files

#### Scenario: Correctness without a GPU backend
- GIVEN a build with no CUDA, OpenCL, or Metal toolkit present
- WHEN a linear `*STATIC` step runs
- THEN the `ComputeBackend` SHALL use the default CPU path
- AND the linear solve SHALL complete and produce correct results without error

### Requirement: Linear static end-to-end pipeline
CalculixPP SHALL execute the complete linear-static pipeline for the supported Phase-1
deck subset, from input parsing to results output, as a single continuous flow. The
supported subset SHALL include `*NODE`, `*ELEMENT` with `C3D4` and `C3D10` tetrahedra,
`*MATERIAL` with `*ELASTIC`, `*SOLID SECTION`, `*BOUNDARY`, `*CLOAD`, `*DLOAD` (pressure),
and `*STATIC`. Element stiffness and load integration SHALL be dispatched through the
`ComputeBackend` (CPU default), and the global solve SHALL use the NumPP sparse direct
solver. Results SHALL be accessible both as written files and through the Python
bindings. (ref: src/CalculiX.c, ref: src/linstatic.c, ref: src/e_c3d.f, ref: src/frd.c)

#### Scenario: Parse-assemble-solve-recover-write walkthrough
- GIVEN an input deck using only the Phase-1 subset (`*NODE`, `*ELEMENT` C3D4/C3D10, `*MATERIAL`/`*ELASTIC`, `*SOLID SECTION`, `*BOUNDARY`, `*CLOAD`, `*DLOAD` pressure, `*STATIC`)
- WHEN the job is run
- THEN CalculixPP SHALL parse the deck into an in-memory model of nodes, tetrahedral elements, material, section, boundary conditions, and loads
- AND assemble the global stiffness `K` and load vector `f` via the `ComputeBackend` CPU path, applying `*BOUNDARY` constraints, `*CLOAD` concentrated loads, and `*DLOAD` pressure loads
- AND solve `K u = f` with the NumPP sparse direct solver
- AND recover element stresses and strains at integration/nodal points from `u`
- AND write nodal displacements and element stresses/strains to `jobname.frd` and requested tabular results to `jobname.dat`

#### Scenario: Results reachable from Python
- GIVEN a completed linear-static run driven through the Python bindings
- WHEN the caller queries the solution
- THEN nodal displacements and element stresses/strains SHALL be retrievable through the Python API in addition to the written `.frd`/`.dat` files

### Requirement: Reference-result agreement
The linear-static results produced by CalculixPP SHALL match the reference CalculiX
results for the equivalent `test/` deck within a documented numerical tolerance. The
acceptance metric SHALL be the relative L2 norm of the nodal displacement field
(`||u_cxx − u_ref||_2 / ||u_ref||_2`) below a documented threshold, with an analogous
tolerance documented for recovered stresses. This agreement SHALL hold on the CPU
backend, which is the correctness oracle. (ref: reference CalculiX `test/` corpus)

#### Scenario: Displacement L2 within tolerance
- GIVEN a Phase-1 `test/` deck and its reference CalculiX result
- WHEN CalculixPP solves the same deck on the CPU backend
- THEN the relative L2 norm of the displacement difference SHALL be at or below the documented tolerance
- AND recovered stresses SHALL agree with the reference within the documented stress tolerance

### Requirement: Nonlinear static analysis
CalculixPP SHALL solve a `*STATIC` step incrementally with Newton-Raphson iterations
when geometric nonlinearity (`*STEP, NLGEOM`), nonlinear material, or contact is
present, updating the tangent stiffness and residual and converging each increment
before advancing step time. Tangent assembly and each increment's solve SHALL route
through the `ComputeBackend` (CPU default) and the NumPP solver. (ref: src/nonlingeo.c)

#### Scenario: Newton-Raphson iteration
- GIVEN a nonlinear `*STATIC` step (`*STEP, NLGEOM`, nonlinear material, or contact)
- WHEN an increment is solved
- THEN CalculixPP SHALL iterate, reassembling the tangent stiffness and residual each iteration
- AND continue until the residual satisfies the convergence controls before advancing time

### Requirement: Incrementation control
A `*STATIC` step SHALL accept an initial time increment, total step time, and
minimum/maximum increments. With automatic incrementation, the increment size SHALL be
reduced on non-convergence and may grow on easy convergence; with `DIRECT` the fixed
increment SHALL be used. Convergence and cutback behavior SHALL be governed by
`*CONTROLS`. (ref: src/nonlingeo.c, ref: src/controlss.f)

#### Scenario: Cutback on divergence
- GIVEN automatic incrementation
- WHEN an increment fails to converge within the allowed iterations
- THEN the increment SHALL be restarted with a smaller time increment
- AND the step SHALL abort if the increment falls below the configured minimum

#### Scenario: Fixed increment with DIRECT
- GIVEN a `*STATIC, DIRECT` step with a specified increment size
- WHEN the step runs
- THEN CalculixPP SHALL use the fixed increment without automatic resizing

### Requirement: Perturbation and buckling preload
CalculixPP SHALL support `*STATIC, PERTURBATION` linear perturbation about a preloaded
base state, providing the stress state used by subsequent perturbation steps and
including its stress stiffening. (ref: src/CalculiX.c, ref: src/nonlingeo.c)

#### Scenario: Preloaded perturbation
- GIVEN a preceding loaded state and a `*STATIC, PERTURBATION` step
- WHEN the step runs
- THEN the response SHALL be linearized about the preloaded base state, including its stress stiffening
