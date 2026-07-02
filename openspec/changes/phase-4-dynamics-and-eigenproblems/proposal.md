## Why

Phases 1–3 solve quasi-static and steady problems. Real structures vibrate, resonate, and
buckle. Phase 4 adds the frequency domain and the time-dependent inertial response: natural
frequencies and mode shapes, linear buckling factors, complex/damped modes, direct and
modal-superposition dynamics, and superelement reduction. The physics is already in the
baseline; the one genuinely-new piece of machinery is the modal/eigensolution engine that
every one of these procedures consumes.

## What Changes

- **Eigenproblems** — `*FREQUENCY` (natural frequencies + mode shapes), cyclic-symmetric
  frequencies, `*BUCKLE` (linear buckling factors about a preloaded base state), and
  `*COMPLEX FREQUENCY` (damped / friction-induced complex modes).
- **Dynamics** — `*DYNAMIC` direct time integration (Newmark / HHT-α), `*MODAL DYNAMIC`
  superposition, `*STEADY STATE DYNAMICS` harmonic response, Rayleigh/modal damping, base
  motion, and the `*GREEN` Green-function step.
- **Substructure / superelement** — `*SUBSTRUCTURE GENERATE` Craig-Bampton reduction using
  fixed-interface normal modes + constraint modes from the eigensolution engine.
- **Extends** `linear-algebra-and-solvers` (NumPP eigensolvers), `constraints` (cyclic
  symmetry), and reuses `nonlinear-solution-control` (prestressed / nonlinear direct
  dynamics) and the Phase-2 incrementation engine (time stepping).
- **New capability** — `eigensolution`: eigenpair extraction, mass-normalized modal basis,
  modal superposition projection, and participation/effective-mass, layered on the raw NumPP
  eigensolve.
- No breaking changes. GPU backends remain optional; advanced physics (Phase 5) deferred.

## Capabilities

### New Capabilities

- `eigensolution`: the reusable modal engine — extract the generalized symmetric eigenpairs
  (via `linear-algebra-and-solvers`), mass-normalize them into a modal basis, compute
  participation factors / effective mass, and project systems and loads onto modal
  coordinates for superposition. Peer to `linear-algebra-and-solvers` and
  `nonlinear-solution-control` (engine, not procedure); consumed by frequency, buckling,
  complex-frequency, modal/steady-state dynamics, and substructure reduction so extraction
  happens once.

### Modified Capabilities

None. `modal-and-buckling-analysis`, `dynamic-analysis`, and `substructure-generation`
already exist as baseline capabilities; this change implements them. `tasks.md` cites the
exact baseline requirements. Only `eigensolution` is added as a delta.

## Impact

- **New code:** the eigensolution engine (subspace/Lanczos with spectral shift on top of
  NumPP, mass normalization, participation, modal projection), the Newmark/HHT direct
  integrator, buckling prestress path, complex-eigen path, cyclic-symmetry sector
  eigenproblem, and Craig-Bampton reduction.
- **Extended code:** parser coverage for the dynamics/eigen card set; Python bindings +
  API-parity check; the regression corpus grows with frequency, buckling, complex-freq,
  direct/modal/steady-state dynamics, and substructure decks.
- **Dependencies:** unchanged (NumPP, CyberCadKernel, pybind11; no GPU required).
- **Depends on:** `phase-1-foundation` (assembly/solve spine) and
  `phase-2-nonlinear-statics-and-materials` (nonlinear driver, incrementation, prestress).
