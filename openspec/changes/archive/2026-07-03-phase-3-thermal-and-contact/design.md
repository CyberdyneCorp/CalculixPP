## Context

Phases 1–2 deliver the assembly/solve spine and the nonlinear driver. Phase 3 introduces
interaction: heat flows between and within bodies, and surfaces touch. Both reuse the
Phase-2 `nonlinear-solution-control` driver — transient conduction, coupled
thermomechanics, and contact are all nonlinear increments with different tangent/residual
callbacks. The one new piece of infrastructure is the geometric contact-search engine,
which mechanical and thermal contact share.

## Goals / Non-Goals

**Goals:**

- Steady + transient conduction, film, cavity radiation, coupled temperature-displacement.
- Node-to-surface and surface-to-surface (mortar) contact with friction and surface behavior.
- Thermal contact (gap conductance, gap heat generation) sharing one contact search.
- Element/contact-pair birth-death (`*MODEL CHANGE`).
- The `contact-search` engine as a reusable capability.

**Non-Goals:**

- Dynamics/eigen (Phase 4); CFD conjugate heat transfer beyond 1-D network film coupling
  (Phase 5); crack propagation (Phase 5).
- GPU contact/thermal kernels (optional, later where hot).
- Wear, adhesion, or cohesive-zone contact (out of scope).

## Decisions

- **contact-search is an engine capability, not part of `contact`.** Search, projection,
  mortar integration, and gap/pressure evaluation live in `contact-search`; the `contact`
  capability specs the user-facing modeling (pairs, surface behavior, friction, output)
  that consumes it. Rationale: the same engine serves mechanical and thermal contact and
  keeps the geometric complexity isolated (systems-band cognitive complexity). Mirrors
  `linear-algebra-and-solvers` / `nonlinear-solution-control`.

- **Mortar dual-Lagrange default for surface-to-surface.** Non-matching meshes use dual
  Lagrange multipliers (condensable, keeps the system SPD-friendly), standard Lagrange
  selectable. Rationale: matches CalculiX mortar behavior and the Phase-1 symmetric solver
  path. Trade-off: dual basis math is intricate; verified against reference mortar decks.

- **Normal contact via penalty with optional augmentation.** Start with a penalty
  (pressure-overclosure) enforcement reusing the surface-behavior law; add augmented
  Lagrange for pressure accuracy where penalty chatters. Alternative: pure Lagrange —
  rejected initially (indefinite systems break the SPD direct path).

- **Coupled thermomechanics: monolithic tangent, staggered fallback.** Assemble the
  coupled temperature-displacement tangent as one system by default (robust strong
  coupling); allow a staggered scheme for weakly coupled, large problems. Rationale:
  monolithic converges best for contact-thermal; staggered scales. Both drive off
  `nonlinear-solution-control`.

- **Radiation view factors computed once per configuration.** Cavity view factors are
  computed geometrically (with symmetry/decimation) and the radiation exchange assembled
  as a dense/sparse coupling solved through NumPP; recomputed only when geometry changes.
  Rationale: view-factor cost is O(n²); avoid per-increment recompute.

- **Transient time integration via the incrementation engine.** Backward-Euler thermal
  transient reuses the Phase-2 automatic incrementation/cutback. Rationale: one time-stepping
  engine across procedures.

## Risks / Trade-offs

- **Contact chattering / non-convergence** → augmented Lagrange, adaptive penalty, and
  cutback from `nonlinear-solution-control`; per-deck limits.
- **Dual-mortar correctness on non-matching meshes** → validate against reference mortar
  decks before enabling by default; standard Lagrange fallback.
- **View-factor cost/accuracy** → symmetry + clustering; document the accuracy knob per deck.
- **Coupling convergence (monolithic vs staggered)** → default monolithic; flag staggered
  as opt-in with a documented convergence caveat.

## Migration Plan

Additive on Phases 1–2. Order: contact-search engine (unit-tested on analytic
projections) → mechanical contact → heat transfer (conduction→film→radiation) → coupled
thermomechanics → model change. Each lands behind the growing regression corpus; new paths
feature-flagged for rollback.

## Open Questions

- Mortar: is dual-Lagrange parity with CalculiX achievable for all reference decks, or is
  standard Lagrange needed for some?
- Radiation solver: dense view-factor matrix vs. iterative — threshold by cavity size?
- Default coupling scheme per problem class (contact-thermal → monolithic; large thermal-only
  → staggered?).
- Phase-3 corpus: which reference decks (steady/transient conduction, film, cavity radiation,
  n2f contact, mortar contact, friction, coupled temp-disp, model change)?
