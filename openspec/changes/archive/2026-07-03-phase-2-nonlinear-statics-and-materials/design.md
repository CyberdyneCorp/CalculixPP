## Context

Phase 1 delivers assembly, a NumPP SPD direct solve, and the linear-static pipeline on
tets. Phase 2 adds the machinery that makes the response nonlinear and the model library
broad enough for real decks. The behavior is fixed by the baseline specs; the design work
is the internal architecture that keeps element math, constitutive integration, and the
solution driver decoupled and testable, and that reuses the Phase-1 ComputeBackend/NumPP
spine without a rewrite.

## Goals / Non-Goals

**Goals:**

- A reusable nonlinear solution driver (Newton-Raphson) shared by later phases.
- Return-mapping plasticity and finite-strain hyperelastic/creep kernels behind the
  Phase-1 material-model interface.
- Hex/wedge element families sharing one isoparametric kernel structure with the tets.
- Constraint elimination (equations, MPC, coupling, tie) folded into assembly via NumPP.
- Automatic incrementation with cutback driven by convergence controls.
- Python bindings + API-parity check extended to every new public surface.

**Non-Goals:**

- Contact (Phase 3), thermal/coupled (Phase 3), dynamics/eigen (Phase 4).
- GPU kernels for the new element/material paths (optional, later where hot).
- Arc-length/Riks solver and material damage *evolution* (initiation only this phase).

## Decisions

- **Nonlinear driver separate from procedure — now a specced capability.** A
  `NonlinearSolver` owns the Newton-Raphson loop (tangent assembly → residual → NumPP
  solve → update → convergence test) and the incrementation engine; `static-analysis`
  supplies the residual/tangent callbacks. This change promotes it to the new
  `nonlinear-solution-control` capability (peer to `linear-algebra-and-solvers`) so the
  driver is specified once and reused. Rationale: Phases 3–4 (heat, dynamics) reuse the
  same driver with different callbacks. Alternative: a static-only loop — rejected, it
  would be re-implemented thrice.

- **Consistent tangent, not just secant.** Element and material kernels return the
  algorithmically consistent tangent so Newton keeps quadratic convergence. Rationale:
  cutback cost dominates if convergence degrades. Trade-off: more kernel math; mitigated
  by finite-difference tangent checks in unit tests.

- **Return-mapping plasticity with per-point state.** Rate-independent plasticity uses a
  radial-return integrator storing state per integration point via the Phase-1 `*DEPVAR`
  mechanism. Rationale: standard, robust, matches CalculiX. Alternative: sub-stepping
  explicit — rejected (accuracy/stability).

- **One isoparametric element kernel, many topologies.** Shape-function/Jacobian
  machinery is topology-parameterized; C3D8/C3D20/C3D6/C3D15 plug into the same
  integration driver as C3D4/C3D10, with selectable reduced integration + hourglass
  control. Rationale: avoids per-element duplication and keeps cognitive complexity in
  the systems band.

- **Constraints eliminated at assembly.** Linear constraints (`*EQUATION`, MPC, coupling,
  tie) build a transformation applied to the assembled system before the NumPP solve,
  rather than Lagrange multipliers, for the SPD-preserving path. Rationale: keeps the
  Phase-1 symmetric direct solver usable. Trade-off: mortar tie needs care; over-constraint
  detection required (already a baseline requirement).

- **Incompressible hyperelastic.** Near-incompressible response uses a mixed
  displacement/pressure (u/p) formulation on the supporting element types. Flagged as a
  risk (below) since it touches the assembly signature.

- **Amplitude engine.** A time-value evaluator (step, tabular, periodic) shared by loads
  and BCs, sampled by the incrementation engine. Rationale: one place for time variation.

## Risks / Trade-offs

- **Newton non-convergence on stiff plasticity** → automatic cutback + line search;
  per-deck iteration/cutback limits from `*CONTROLS`; quarantine pathological decks.
- **Mixed u/p changes the element interface** → introduce the u/p path behind the same
  kernel abstraction; keep pure-displacement elements unaffected; gate behind material
  incompressibility.
- **Constraint elimination vs. SPD assumption** → validate that eliminated systems stay
  symmetric; fall back to the general solver path (Phase-1 unavailable-solver policy) if not.
- **Regression tolerance drift on path-dependent results** → per-deck tolerance manifest
  (from Phase 1), stress compared at integration points, displacement L2 primary.

## Migration Plan

Additive on top of `phase-1-foundation`; no data migration. Land the nonlinear driver
first (validated on a linear deck → identical to the Phase-1 direct solve in one
increment), then materials, elements, loads, constraints, each behind the growing
regression corpus. Rollback is per-workstream (feature-flag new element/material paths).

## Open Questions

- Line-search: needed for Phase-2 material set, or does cutback suffice initially?
- Reduced integration + hourglass control: which stabilization (assumed-strain vs.
  physical hourglass) to match CalculiX results within tolerance?
- Which reference `test/` decks form the Phase-2 corpus (one per: plasticity, hyperelastic,
  creep, MPC/equation, coupling, tie, hex/wedge, body load)?
- Mortar `*TIE` for non-matching meshes — reuse contact-mortar infrastructure (Phase 3) or
  a standalone tie projection now?
