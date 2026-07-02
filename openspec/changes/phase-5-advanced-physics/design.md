## Context

Phases 1–4 deliver the structural core (statics, thermal/contact, dynamics) and its engine
layer. Phase 5 adds the remaining physics and the adaptivity/reuse capabilities. Three of
these procedures (CFD, electromagnetic-thermal, network-thermal) share the same need: solve
several fields and exchange coupling terms until the coupled system converges. That shared
need becomes the `field-coupling` engine, completing the engine layer (linear solve →
nonlinear driver → contact search → eigensolution → field coupling). The other Phase-5
capabilities reuse the existing spine and CyberCadKernel.

## Goals / Non-Goals

**Goals:**

- 3-D CFD + 1-D fluid networks and their coupling to heat transfer.
- Electromagnetics (magnetostatics, eddy currents, electric conduction) + Joule coupling.
- Crack propagation (SIF, growth, remeshing), design sensitivity/optimization, adaptive
  mesh refinement, submodeling, and high-cycle fatigue.
- The `field-coupling` engine as a reusable capability.

**Non-Goals:**

- Turbulence modeling beyond what the reference solver provides (match reference scope).
- New GPU kernels for CFD/EM (optional, later).
- Extracting adjoint-sensitivity or remeshing as separate engine capabilities this phase
  (they live inside their owning capabilities to bound scope).

## Decisions

- **`field-coupling` is the single new engine.** Segregated iteration, coupling-term
  exchange, pressure-velocity coupling, and coupled convergence live here; the physics
  procedures (`cfd-and-network-analysis`, `electromagnetic-analysis`,
  `heat-transfer-analysis`) supply the field solves and coupling terms. Rationale: one
  place for multi-physics orchestration; mirrors the earlier engine capabilities.

- **Segregated, not monolithic, multi-physics.** Fields are solved in sequence with
  under-relaxed coupling-term exchange. Rationale: matches the reference CFD/EM approach and
  keeps each field on its existing solver; monolithic multi-physics would rebuild every
  assembly path. Trade-off: outer-loop convergence needs relaxation control (specced).

- **SIMPLE-class pressure-velocity coupling for CFD.** Rationale: standard, robust for
  incompressible/compressible internal flows; solved through NumPP.

- **Adjoint sensitivity inside `design-optimization`.** The adjoint solve reuses the
  factorized system from the primal analysis (NumPP), producing gradients for
  objectives/constraints. Rationale: tightly bound to the optimization procedure; not
  reused elsewhere enough to justify a separate engine this phase.

- **Crack + refinement remeshing via CyberCadKernel.** Crack-front advance and `*REFINE MESH`
  both use CyberCadKernel tet remeshing (reusing the `mesh-processing` quality metrics); the
  solve–refine–resolve loop reuses the analysis driver. Rationale: one remeshing path.

- **HCF and submodeling reuse prior results.** HCF post-processes Phase-4 dynamic/modal
  results; submodeling interpolates a global run's boundary results onto the local model
  (reusing the `mesh-processing`/results interpolation). Rationale: no new solve machinery.

## Risks / Trade-offs

- **CFD stability / convergence** → conservative under-relaxation defaults, coupled residual
  stall detection, per-deck tuning; keep turbulence scope to the reference.
- **EM gauge / air-domain conditioning** → follow the reference formulation; validate on EM
  reference decks; iterative NumPP solve with appropriate preconditioner.
- **Crack remeshing robustness** → reuse `mesh-processing` sliver removal + quality gates;
  fall back / report when a remesh degrades below quality thresholds.
- **Adjoint gradient correctness** → finite-difference gradient checks in the test suite.
- **Coupling false convergence** → outer tolerance separate from inner; stall detection.

## Migration Plan

Additive on Phases 1–4. Order: field-coupling engine (validated on a simple two-field
manufactured coupling) → CFD/networks → electromagnetics + Joule → crack → sensitivity/
optimization → adaptive refinement → submodeling → HCF. Each behind the growing regression
corpus; new paths feature-flagged. This change completes the port; on archive the baseline
gains `field-coupling` and every recognized card has an implemented, validated home.

## Open Questions

- CFD turbulence scope: which models (if any) does the reference corpus exercise?
- Default coupling scheme and relaxation per problem class (CFD vs EM-thermal vs network).
- Phase-5 corpus: representative reference decks per procedure (CFD, network, EM, crack,
  optimization, refinement, submodel, HCF).
- Optimization: which optimizer(s) drive `*FEASIBLE DIRECTION` and robust design to match
  the reference within tolerance.
