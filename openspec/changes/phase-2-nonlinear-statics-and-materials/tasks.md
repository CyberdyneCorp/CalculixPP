# Phase 2 — Nonlinear Statics & Materials

Implementation against the front-loaded baseline plus the one new capability
(`nonlinear-solution-control`) added by this change. Each task cites the capability it
satisfies. Depends on `phase-1-foundation`.

## 1. Nonlinear solution driver (spec: nonlinear-solution-control [NEW]; static-analysis — "Nonlinear static analysis", "Incrementation control")

- [x] 1.1 `NonlinearSolver` Newton-Raphson loop: tangent assembly → residual → solve → update (`numerics::solve_nonlinear_static`; tangent = `assemble_linear_static` reduced system, residual via extracted `fem::internal_force`, solve via `solve_reduced`/`resolve_solver_kind`)
- [x] 1.2 Residual/convergence controls (force + displacement) configurable via `*CONTROLS` (`NonlinearControls`, documented defaults; parser reads `PARAMETERS=FIELD` and `TIME INCREMENTATION`)
- [x] 1.3 Automatic incrementation engine with cutback on divergence; `DIRECT` fixed increment; `*TIME POINTS` (`Incrementation`, `run_increments`/`clamp_increment`; grows after easy convergence, halves on failure, aborts below min)
- [x] 1.4 Line search (optional, behind a flag) for difficult increments (`NonlinearOptions::line_search`, OFF by default; scales the Newton update to reduce the residual)
- [x] 1.5 Validation: a linear deck solved in one increment reproduces the Phase-1 direct solve exactly (single-tet + beam10p reproduce `solve_linear_static` to rel-L2 < 1e-10; `tests/test_nonlinear.cpp`, `python/tests/test_regression.py`)
- [~] 1.6 `*STATIC, PERTURBATION` about a preloaded base state (stress stiffening) — DEFERRED: the geometric stress-stiffening K_geo needs geometric nonlinearity (NLGEOM) that does not exist yet; only linear elasticity is implemented, so a faithful perturbation-about-preload solve cannot be built without faking K_geo. Revisit with workstream 3/4 (finite-strain elements + plasticity).

## 2. Amplitude & load breadth (spec: loads-and-boundary-conditions)

- [x] 2.1 Amplitude engine: step, tabular, and periodic `*AMPLITUDE`, sampled per increment (`core/amplitude.hpp` `Amplitude::value_at`; `Model::amplitude_factor`; `*CLOAD`/`*DLOAD`/`*BOUNDARY` carry `AMPLITUDE=`; driver samples per-increment via `fem::external_load_vector(model, lambda)` and amplitude-scaled prescribed BCs; `tests/test_loads.cpp`)
- [x] 2.2 Body loads: gravity and centrifugal `*DLOAD`; `*DSLOAD` (`BodyLoad` GRAV/CENTRIF integrated as `rho*N` over element volume via the shape/Gauss rule; `*DSLOAD` reuses the pressure-face machinery; `fem::external_load_vector` extended; gravity global-equilibrium + centrifugal outward-radial tests in `tests/test_loads.cpp`)
- [~] 2.3 Step/restart property changes: `*CHANGE MATERIAL`, `*CHANGE PLASTIC`, `*CHANGE SOLID SECTION`. `*CHANGE SOLID SECTION` (MATERIAL=/ELSET=) is DONE — parsed and applied by appending a `SolidSection` that wins per element (`element_elastic` is last-writer); re-binds material within the single step. `*CHANGE MATERIAL, NAME=` is parsed (re-opens the named material, rejects unknown names) and `*CHANGE PLASTIC` data is accepted, but PARTIAL/DEFERRED: it only ever changes plastic hardening data, and plasticity does not exist yet (workstream 4), so the change is a no-op. The *cross-step / restart* aspect (redefining properties at a step boundary) also needs multi-step step handling, which the single-step model lacks. Parser + within-step apply for solid-section rebind implemented and tested (`tests/test_step_changes.cpp`).
- [~] 2.4 Load accumulation semantics (`OP=MOD` / `OP=NEW`). Parsed and applied *within the single step*: `OP=` on `*CLOAD`/`*DLOAD`/`*BOUNDARY` is validated (MOD default / NEW), and `OP=NEW` resets the accumulated loads of that type (concentrated / distributed+body / SPC), honoring the "first card of the type only" rule via per-type once-per-step flags. DEFERRED: true *across-steps* semantics — OP=MOD keeping prior-step loads while adding same-step same-DOF forces, and OP=NEW removing prior-step loads — require multi-step step handling that does not exist yet (the model is single-step; there is no step boundary to accumulate across). Within-step reset implemented and tested (`tests/test_step_changes.cpp`).

## 3. Element library expansion (spec: element-sections, mesh-and-model)

- [ ] 3.1 Shared isoparametric kernel parameterized by topology (shape fns, Jacobian, integration)
- [ ] 3.2 Hex/wedge families: `C3D8`, `C3D8R`, `C3D20`, `C3D20R`, `C3D6`, `C3D15`
- [ ] 3.3 Reduced integration + hourglass/stabilization control
- [ ] 3.4 Shell / beam / membrane sections (`*SHELL SECTION`, `*BEAM SECTION`, `*MEMBRANE SECTION`)
- [ ] 3.5 Discrete/connector elements: `*MASS`, `*SPRING`, `*DASHPOT`
- [ ] 3.6 `*USER ELEMENT` / `*USER SECTION` extensibility; `*NODAL THICKNESS` / `*NORMAL`; `*DISTRIBUTION`

## 4. Constitutive kernels (spec: material-models)

- [ ] 4.1 Return-mapping rate-independent plasticity: `*PLASTIC` isotropic + kinematic, `*CYCLIC HARDENING`
- [ ] 4.2 Consistent tangent for plasticity; finite-difference tangent verification test
- [ ] 4.3 Hyperelasticity / foam: `*HYPERELASTIC`, `*HYPERFOAM` (with near-incompressible u/p path)
- [ ] 4.4 Creep / visco: `*CREEP`, `*VISCO`, `*VALUES AT INFINITY`
- [ ] 4.5 `*DEFORMATION PLASTICITY`, `*MOHR COULOMB` + `*MOHR COULOMB HARDENING`, `*DAMAGE INITIATION`
- [ ] 4.6 User material: `*USER MATERIAL` C++20 interface, `*DEPVAR` state, `*RATEDEPENDENT` scaling

## 5. Constraints (spec: constraints)

- [ ] 5.1 `*EQUATION` linear constraints eliminated at assembly (SPD-preserving transform)
- [ ] 5.2 `*MPC` PLANE / STRAIGHT / BEAM (and user-defined hook)
- [ ] 5.3 `*RIGID BODY`
- [ ] 5.4 `*COUPLING` kinematic + distributing / `*DISTRIBUTING COUPLING`
- [ ] 5.5 Mortar `*TIE` for matching + non-matching surfaces; over-constraint detection

## 6. Parser & bindings & validation

- [ ] 6.1 Extend the deck parser to the full Phase-2 card set
- [ ] 6.2 Extend Python bindings to all new public surfaces; extend the API-parity coverage test
- [ ] 6.3 Curate the Phase-2 reference corpus (plasticity, hyperelastic, creep, MPC/equation, coupling, tie, hex/wedge, body load) + per-deck tolerances
- [ ] 6.4 pytest regression: each Phase-2 deck matches reference CalculiX within tolerance
- [ ] 6.5 Update README/docs with a nonlinear example
