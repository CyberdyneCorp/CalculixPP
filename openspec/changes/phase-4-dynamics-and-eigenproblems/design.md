## Context

Phases 1–3 deliver the assembly/solve spine, the nonlinear driver, and interaction
(contact/thermal). Phase 4 adds inertia and the frequency domain. Every Phase-4
procedure — frequency, buckling, complex frequency, modal/steady-state dynamics, and
Craig-Bampton substructure reduction — needs the same underlying operation: extract a
modal basis and (for the dynamics procedures) superpose in modal coordinates. That shared
operation becomes the `eigensolution` engine, layered on the raw NumPP eigensolve that
`linear-algebra-and-solvers` already specs.

## Goals / Non-Goals

**Goals:**

- `*FREQUENCY`, cyclic-symmetric frequency, `*BUCKLE`, `*COMPLEX FREQUENCY`.
- `*DYNAMIC` direct integration, `*MODAL DYNAMIC`, `*STEADY STATE DYNAMICS`, damping, base
  motion, `*GREEN`.
- `*SUBSTRUCTURE GENERATE` Craig-Bampton reduction (fixed-interface + constraint modes).
- The `eigensolution` engine as a reusable capability.

**Non-Goals:**

- CFD/acoustic fluid dynamics, electromagnetics, crack, optimization (Phase 5).
- GPU eigen/integration kernels (optional, later).
- Nonlinear frequency-domain (harmonic balance) — out of scope.

## Decisions

- **`eigensolution` is an engine on top of the raw eigensolve.** `linear-algebra-and-solvers`
  provides the generalized symmetric eigenpairs from NumPP; `eigensolution` adds mass
  normalization, spectral shift/rigid-body handling, participation/effective mass, modal
  projection, and the complex path — the reusable modal machinery. Rationale: avoids
  duplicating the raw solver while giving the procedures one modal API. Mirrors
  `nonlinear-solution-control` over the linear solve.

- **Subspace iteration / Lanczos with shift-invert.** Default extraction is shift-invert
  Lanczos on a NumPP factorization, robust for rigid-body and buckling (non-positive)
  eigenvalues. Rationale: matches CalculiX/ARPACK behavior and handles the shifted spectrum.
  Alternative: direct dense eigensolve — rejected for size.

- **Newmark / HHT-α for direct integration.** Direct `*DYNAMIC` uses the unconditionally
  stable HHT-α scheme with controllable numerical damping, time-stepped by the Phase-2
  incrementation engine; nonlinear direct dynamics (contact, plasticity) drive off
  `nonlinear-solution-control`. Rationale: one time-stepping engine, one Newton driver.

- **Modal superposition for linear transient/harmonic.** Linear `*MODAL DYNAMIC` and
  `*STEADY STATE DYNAMICS` integrate in modal coordinates (cheap, decoupled); nonlinear or
  contact dynamics stay in physical coordinates via direct integration. Rationale: cost.

- **Buckling is a two-step prestress + eigenproblem.** Compute the base stress state
  (reuse the static prestress path), form the geometric stiffness, and solve
  `(K + λ K_geo) x = 0`. Rationale: standard linear buckling; reuses Phase-1/2 assembly.

- **Cyclic symmetry via a sector complex eigenproblem.** Extract per-nodal-diameter
  eigenpairs on a single sector with complex cyclic constraints (reuse `constraints`
  cyclic-symmetry). Rationale: avoids meshing the full wheel.

- **Craig-Bampton reduction reuses fixed-interface modes.** `substructure-generation`
  requests fixed-interface normal modes from `eigensolution` and combines them with
  constraint modes over the retained DOFs to form the reduced matrices. Rationale: no
  separate eigen path for superelements.

## Risks / Trade-offs

- **Missed / starved modes** → Sturm-sequence check on the shifted factorization to verify
  the count in an interval; widen the subspace on shortfall.
- **Large eigen cost** → shift-invert factorization reused across modes; block Lanczos;
  document mode-count/cost per deck.
- **Complex-eigen robustness** → validate against reference complex-frequency decks; fall
  back to a larger Krylov subspace.
- **HHT-α numerical damping calibration** → default α from CalculiX; expose the knob; verify
  energy behavior on a reference transient.

## Migration Plan

Additive on Phases 1–3. Order: eigensolution engine (validated on an analytic
free-vibration eigenproblem) → frequency → buckling → complex frequency → direct dynamics →
modal/steady-state dynamics → substructure → cyclic symmetry. Each behind the growing
regression corpus; new paths feature-flagged.

## Open Questions

- Eigensolver algorithm parity with CalculiX/ARPACK: Lanczos vs. subspace iteration for
  which problem sizes?
- Default damping model (Rayleigh coefficients vs. modal damping ratios) and its parsing.
- Phase-4 corpus: which reference decks (frequency, cyclic frequency, buckle, complex freq,
  direct dynamic, modal dynamic, steady-state, substructure)?
- Steady-state dynamics output: which frequency sweep / harmonic result set to compare.
