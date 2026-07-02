## Why

Phases 1–2 solve isolated bodies. Real assemblies exchange heat and touch. Phase 3 adds
the first capabilities where separate surfaces and fields interact: thermal conduction and
its coupling to the mechanical response, and contact between surfaces — plus staged
activation/deactivation of parts (`*MODEL CHANGE`). This is where the solver becomes usable
for assemblies, thermal stress, and staged/assembly analyses. The physics is already in the
baseline; the one genuinely-new piece of machinery is the geometric contact-search engine.

## What Changes

- **Heat transfer** — steady and transient conduction, convective film (`*FILM`), cavity
  radiation with view factors (`*RADIATE`), and coupled temperature-displacement
  (`*COUPLED TEMPERATURE-DISPLACEMENT`), reusing the Phase-2 nonlinear driver for the
  transient/coupled increments.
- **Contact** — node-to-surface and surface-to-surface (mortar) formulations, surface
  behavior (hard / exponential / linear / tied pressure-overclosure), friction (Coulomb,
  stick↔slip), thermal contact (`*GAP CONDUCTANCE`, `*GAP HEAT GENERATION`), and contact
  output.
- **Model change** — element and contact-pair birth/death between steps (excavation,
  additive, staged assembly).
- **Extends** `material-models` (conductivity/specific-heat consumed), `results-output`
  (thermal fields `NT`/`HFL`), and the Phase-2 `nonlinear-solution-control` driver
  (transient/coupled/contact iterations).
- **New capability** — `contact-search`: the geometric proximity search + node/face
  projection + mortar segment integration + gap/pressure evaluation engine.
- No breaking changes. GPU backends remain optional; dynamics/eigen (Phase 4) and advanced
  physics (Phase 5) stay deferred.

## Capabilities

### New Capabilities

- `contact-search`: the reusable geometric contact engine — spatial proximity search,
  node-to-surface projection, surface-to-surface mortar segment integration, and
  gap/pressure-overclosure + stick/slip evaluation — shared by mechanical and thermal
  contact. Peer to `linear-algebra-and-solvers` and `nonlinear-solution-control`
  (engine, not procedure); the `contact` capability specs the user-facing modeling that
  consumes it.

### Modified Capabilities

None. `heat-transfer-analysis`, `contact` (incl. thermal-contact requirements), and
`model-change` already exist as baseline capabilities; this change implements them.
`tasks.md` cites the exact baseline requirements. Only `contact-search` is added as a delta.

## Impact

- **New code:** contact-search engine (spatial acceleration, projection, mortar
  integration), thermal element kernels, radiation view-factor computation, the
  thermomechanical coupling scheme, film/radiation surface loads.
- **Extended code:** parser coverage for the thermal/contact card set; Python bindings +
  API-parity check; the regression corpus grows with thermal, contact, and coupled decks.
- **Dependencies:** unchanged (NumPP, CyberCadKernel, pybind11; no GPU required).
- **Depends on:** `phase-1-foundation` (assembly/solve spine) and
  `phase-2-nonlinear-statics-and-materials` (nonlinear driver, incrementation).
