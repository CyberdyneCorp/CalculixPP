# Phase 2 — Nonlinear Statics & Materials

Implementation against the front-loaded baseline plus the one new capability
(`nonlinear-solution-control`) added by this change. Each task cites the capability it
satisfies. Depends on `phase-1-foundation`.

## 1. Nonlinear solution driver (spec: nonlinear-solution-control [NEW]; static-analysis — "Nonlinear static analysis", "Incrementation control")

- [ ] 1.1 `NonlinearSolver` Newton-Raphson loop: tangent assembly → residual → NumPP solve → update
- [ ] 1.2 Residual/convergence controls (force + displacement) configurable via `*CONTROLS`
- [ ] 1.3 Automatic incrementation engine with cutback on divergence; `DIRECT` fixed increment; `*TIME POINTS`
- [ ] 1.4 Line search (optional, behind a flag) for difficult increments
- [ ] 1.5 Validation: a linear deck solved in one increment reproduces the Phase-1 direct solve exactly
- [ ] 1.6 `*STATIC, PERTURBATION` about a preloaded base state (stress stiffening)

## 2. Amplitude & load breadth (spec: loads-and-boundary-conditions)

- [ ] 2.1 Amplitude engine: step, tabular, and periodic `*AMPLITUDE`, sampled per increment
- [ ] 2.2 Body loads: gravity and centrifugal `*DLOAD`; `*DSLOAD`
- [ ] 2.3 Step/restart property changes: `*CHANGE MATERIAL`, `*CHANGE PLASTIC`, `*CHANGE SOLID SECTION`
- [ ] 2.4 Load accumulation semantics (`OP=MOD` / `OP=NEW`) across steps

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
