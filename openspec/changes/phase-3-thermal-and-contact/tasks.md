# Phase 3 — Thermal & Contact

Implements baseline thermal/contact/model-change physics plus the one new capability
(`contact-search`). Each task cites the capability it satisfies. Depends on
`phase-1-foundation` and `phase-2-nonlinear-statics-and-materials`.

## 1. Contact-search engine (spec: contact-search [NEW])

- [ ] 1.1 Spatial proximity search (bounding structure) with configurable search/adjust distance
- [ ] 1.2 Node-to-surface projection: closest point, signed gap, local normal/tangent frame
- [ ] 1.3 Surface-to-surface mortar segment integration (dual + standard Lagrange basis)
- [ ] 1.4 Gap / pressure-overclosure evaluation and tangential stick/slip state
- [ ] 1.5 Contact operator generation assembled via ComputeBackend/NumPP within the Newton iteration
- [ ] 1.6 Expose gap/pressure/area state for thermal contact reuse

## 2. Contact modeling (spec: contact)

- [ ] 2.1 `*CONTACT PAIR` (node-to-surface + surface-to-surface) and `*SURFACE INTERACTION`
- [ ] 2.2 `*SURFACE BEHAVIOR` (hard / exponential / linear / tied pressure-overclosure)
- [ ] 2.3 `*FRICTION` (Coulomb, stick stiffness, stick↔slip)
- [ ] 2.4 Thermal contact: `*GAP CONDUCTANCE`, `*GAP HEAT GENERATION`
- [ ] 2.5 Contact modifiers + output (`*CLEARANCE`, `*CONTACT FILE`/`PRINT`/`OUTPUT`, CSTR)

## 3. Heat transfer (spec: heat-transfer-analysis)

- [ ] 3.1 Steady-state conduction element kernels; `*HEAT TRANSFER, STEADY STATE`
- [ ] 3.2 Transient conduction (backward Euler via the incrementation engine)
- [ ] 3.3 Convective film loads (`*FILM`) and forced-convection network coupling
- [ ] 3.4 Cavity radiation: view-factor computation + radiation exchange assembled via NumPP (`*RADIATE`)

## 4. Coupled thermomechanics (spec: heat-transfer-analysis — coupled)

- [ ] 4.1 `*COUPLED TEMPERATURE-DISPLACEMENT` monolithic tangent through nonlinear-solution-control
- [ ] 4.2 Staggered coupling fallback (opt-in) for large weakly-coupled problems
- [ ] 4.3 Thermal expansion coupling verified against reference thermal-stress decks

## 5. Model change (spec: model-change)

- [ ] 5.1 `*MODEL CHANGE, TYPE=ELEMENT` deactivate / reactivate element sets (strain-free reactivation)
- [ ] 5.2 `*MODEL CHANGE, TYPE=CONTACT PAIR` activate / deactivate contact pairs between steps

## 6. Materials, results, parser, bindings, validation

- [ ] 6.1 Consume `*CONDUCTIVITY` / `*SPECIFIC HEAT` in thermal kernels (spec: material-models)
- [ ] 6.2 Thermal result fields `NT` / `HFL` in `.frd` / `.dat` (spec: results-output)
- [ ] 6.3 Extend deck parser to the Phase-3 card set
- [ ] 6.4 Extend Python bindings + API-parity coverage test
- [ ] 6.5 Phase-3 reference corpus (steady/transient conduction, film, radiation, n2f + mortar contact, friction, coupled temp-disp, model change) + per-deck tolerances
- [ ] 6.6 pytest regression: each Phase-3 deck matches reference CalculiX within tolerance
- [ ] 6.7 Update README/docs with a thermal + contact example
