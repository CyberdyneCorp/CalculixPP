## Why

Phase 1 proved the linear-static pipeline end to end. It does not yet solve the
problems real decks contain: nonlinear response (large deformation, plasticity,
contact-free but path-dependent behavior), the wider element and material libraries, and
the constraints that tie real assemblies together. Phase 2 turns the working skeleton
into a usable nonlinear structural solver. All Phase-2 behavior is already specified in
the front-loaded baseline; this change **implements** that subset — it does not add or
change requirements.

## What Changes

- **Nonlinear static solution** — Newton-Raphson with a consistent tangent, residual
  convergence controls, automatic time incrementation with cutback (`*CONTROLS`,
  `*TIME POINTS`), and `*STATIC, PERTURBATION` about a preloaded base state.
- **Material breadth** — rate-independent plasticity (`*PLASTIC` isotropic/kinematic,
  `*CYCLIC HARDENING`), hyperelasticity/foam (`*HYPERELASTIC`, `*HYPERFOAM`), creep/visco
  (`*CREEP`, `*VISCO`, `*VALUES AT INFINITY`), `*DEFORMATION PLASTICITY`, `*MOHR COULOMB`
  (+`*MOHR COULOMB HARDENING`), `*DAMAGE INITIATION`, and user materials via the C++20
  interface (`*USER MATERIAL`, `*DEPVAR`, `*RATEDEPENDENT`).
- **Element breadth** — hex/wedge families (`C3D8(R)`, `C3D20(R)`, `C3D6`, `C3D15`) beside
  the Phase-1 tets; shell/beam/membrane sections; discrete `*MASS`/`*SPRING`/`*DASHPOT`;
  `*USER ELEMENT`/`*USER SECTION`; `*NODAL THICKNESS`/`*NORMAL`; `*DISTRIBUTION`.
- **Loads breadth** — distributed/body loads (gravity, centrifugal), `*DSLOAD`, tabular
  and periodic `*AMPLITUDE`, and step/restart property changes (`*CHANGE MATERIAL`,
  `*CHANGE PLASTIC`, `*CHANGE SOLID SECTION`).
- **Constraints** — `*EQUATION`, `*MPC` (PLANE/STRAIGHT/BEAM), `*RIGID BODY`,
  `*COUPLING` (kinematic + distributing), and mortar `*TIE`, eliminated/enforced during
  assembly through NumPP.
- No breaking changes. GPU backends remain optional; thermal, dynamics, and advanced
  physics stay deferred to Phases 3–5.

## Capabilities

### New Capabilities

- `nonlinear-solution-control`: the reusable Newton-Raphson driver, convergence controls
  (`*CONTROLS`), automatic incrementation/cutback (`*TIME POINTS`, `DIRECT`), amplitude
  time-stepping, and the consistent-tangent contract. This is the nonlinear analog of the
  existing `linear-algebra-and-solvers` capability — the iterative-solve engine as a
  capability in its own right, previously only implicit inside the procedure specs, and
  reused by Phases 3–4.

### Modified Capabilities

None. The physics breadth implemented here (nonlinear statics, plasticity, element/load
families, constraints) already exists as (Phase-2-tagged) requirements in the baseline
specs — `static-analysis`, `material-models`, `element-sections`,
`loads-and-boundary-conditions`, and `constraints`. Implementing an existing requirement
is not a spec-level change, so those carry no delta; `tasks.md` cites the exact baseline
requirements being satisfied. Only the genuinely-new `nonlinear-solution-control`
capability is added as a delta.

## Impact

- **New code:** nonlinear solution driver (Newton-Raphson, line search, incrementation
  engine), constitutive kernels (return-mapping plasticity, hyperelastic/creep), the
  hex/wedge element families, constraint assembly/elimination, and the amplitude engine.
- **Extended code:** parser coverage for the Phase-2 card set; Python bindings for the
  new surfaces (with the API-parity check extended); the pytest regression corpus grows
  with nonlinear/plastic/MPC reference decks.
- **Dependencies:** unchanged (NumPP, CyberCadKernel, pybind11; no GPU required).
- **Depends on:** `phase-1-foundation` (build, CPU backend, assembly/solve spine).
