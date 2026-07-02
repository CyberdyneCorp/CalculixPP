# Phase 4 — Dynamics & Eigenproblems

Implements baseline dynamics/eigen/substructure physics plus the one new capability
(`eigensolution`). Each task cites the capability it satisfies. Depends on
`phase-1-foundation` and `phase-2-nonlinear-statics-and-materials`.

## 1. Eigensolution engine (spec: eigensolution [NEW])

- [ ] 1.1 Generalized symmetric eigenpair extraction via NumPP (subspace/Lanczos, shift-invert)
- [ ] 1.2 Mass normalization and ascending-order modal basis
- [ ] 1.3 Spectral shift + rigid-body / non-positive eigenvalue handling; Sturm-sequence count check
- [ ] 1.4 Participation factors and modal effective mass
- [ ] 1.5 Modal superposition projection of operators and loads
- [ ] 1.6 Complex / damped eigenproblem path

## 2. Eigenfrequency, buckling, complex frequency (spec: modal-and-buckling-analysis)

- [ ] 2.1 `*FREQUENCY` natural frequencies + mode shapes
- [ ] 2.2 `*BUCKLE` two-step prestress + geometric-stiffness eigenproblem
- [ ] 2.3 `*COMPLEX FREQUENCY` damped / friction-induced modes
- [ ] 2.4 Preceding-frequency-step requirement enforced for modal procedures

## 3. Direct-integration dynamics (spec: dynamic-analysis)

- [ ] 3.1 Newmark / HHT-α integrator time-stepped by the incrementation engine
- [ ] 3.2 Nonlinear direct dynamics through `nonlinear-solution-control` (contact/plastic inertial)
- [ ] 3.3 Numerical-damping (α) control and energy verification

## 4. Modal & steady-state dynamics (spec: dynamic-analysis)

- [ ] 4.1 `*MODAL DYNAMIC` superposition using the eigensolution basis
- [ ] 4.2 `*STEADY STATE DYNAMICS` harmonic response over a frequency sweep
- [ ] 4.3 Damping (Rayleigh + modal) and base motion
- [ ] 4.4 `*GREEN` Green-function step

## 5. Substructure / superelement (spec: substructure-generation)

- [ ] 5.1 Fixed-interface normal modes from eigensolution + constraint modes over retained DOFs
- [ ] 5.2 Craig-Bampton reduced stiffness/mass assembly (`*SUBSTRUCTURE GENERATE`, `*RETAINED NODAL DOFS`)
- [ ] 5.3 Reduced / global matrix export (`*SUBSTRUCTURE MATRIX OUTPUT`, `*MATRIX ASSEMBLE`)

## 6. Cyclic symmetry (spec: constraints)

- [ ] 6.1 Cyclic-symmetry sector complex eigenproblem (`*CYCLIC SYMMETRY MODEL`)
- [ ] 6.2 Nodal-diameter mode selection (`*SELECT CYCLIC SYMMETRY MODES`)

## 7. Parser, bindings, validation

- [ ] 7.1 Extend the deck parser to the Phase-4 card set
- [ ] 7.2 Extend Python bindings + API-parity coverage test
- [ ] 7.3 Phase-4 reference corpus (frequency, cyclic frequency, buckle, complex freq, direct/modal/steady-state dynamics, substructure) + per-deck tolerances
- [ ] 7.4 pytest regression: each Phase-4 deck matches reference CalculiX within tolerance
- [ ] 7.5 Update README/docs with a modal + dynamics example
