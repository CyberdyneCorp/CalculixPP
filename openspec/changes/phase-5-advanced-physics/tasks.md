# Phase 5 — Advanced Physics & Adaptivity

Implements the remaining baseline physics/adaptivity capabilities plus the one new
capability (`field-coupling`). Each task cites the capability it satisfies. Depends on
`phase-1-foundation` through `phase-4-dynamics-and-eigenproblems`. Completes the port.

## 1. Field-coupling engine (spec: field-coupling [NEW])

- [ ] 1.1 Segregated multi-field iteration with coupling-term exchange and outer convergence
- [ ] 1.2 Coupling-term transfer (Joule heat, temperature-dependent properties, film, tractions)
- [ ] 1.3 SIMPLE-class pressure-velocity coupling for Navier-Stokes
- [ ] 1.4 Coupled convergence control + under-relaxation / stall detection
- [ ] 1.5 Conjugate fluid–solid / network–solid interface flux+state continuity

## 2. CFD & 1-D networks (spec: cfd-and-network-analysis)

- [ ] 2.1 3-D CFD (`*CFD`) compressible/incompressible Navier-Stokes via the field-coupling engine
- [ ] 2.2 1-D fluid-network elements (orifice, pipe, valve, labyrinth) + `*FLUID SECTION` types
- [ ] 2.3 Fluid/physical constants (`*FLUID CONSTANTS`, `*SPECIFIC GAS CONSTANT`, `*PHYSICAL CONSTANTS`)
- [ ] 2.4 Face/flow BCs (`*BOUNDARYF`, `*MASS FLOW`, `*TRANSFORMF`); network MPCs
- [ ] 2.5 Network ↔ heat-transfer coupling (forced-convection film)

## 3. Electromagnetics (spec: electromagnetic-analysis)

- [ ] 3.1 `*ELECTROMAGNETICS` magnetostatics / eddy-current / electric-conduction potential solves
- [ ] 3.2 EM material properties (`*MAGNETIC PERMEABILITY`, `*ELECTRICAL CONDUCTIVITY`)
- [ ] 3.3 Joule-heating coupling to the thermal field via field-coupling
- [ ] 3.4 Enveloping air/surrounding domain

## 4. Crack propagation (spec: crack-propagation)

- [ ] 4.1 `*CRACK PROPAGATION` procedure: SIF along the crack front
- [ ] 4.2 Crack length/shape advance with CyberCadKernel remeshing
- [ ] 4.3 Growth-rate law and crack results output

## 5. Design sensitivity & optimization (spec: design-optimization)

- [ ] 5.1 Design variables/responses/objectives/constraints (`*DESIGN VARIABLES`, `*OBJECTIVE`, `*CONSTRAINT`)
- [ ] 5.2 Adjoint sensitivity reusing the primal NumPP factorization; finite-difference gradient checks
- [ ] 5.3 Filtering / geometric tolerances (`*FILTER`, `*GEOMETRIC TOLERANCE`)
- [ ] 5.4 Feasible-direction optimization (`*FEASIBLE DIRECTION`)
- [ ] 5.5 Robust design / random field (`*ROBUST DESIGN`, `*RANDOM FIELD`, `*CORRELATION LENGTH`)

## 6. Adaptive mesh refinement (spec: mesh-refinement)

- [ ] 6.1 `*REFINE MESH` field-driven tet refinement via CyberCadKernel; sliver removal
- [ ] 6.2 Solve–refine–resolve loop (automatic restart on the refined mesh)
- [ ] 6.3 `*INITIAL MESH` / `*USE REFINED MESH`; refined-mesh writeback (quadratic order preserved)

## 7. Submodeling (spec: submodeling)

- [ ] 7.1 `*SUBMODEL` boundary declaration + global-result source
- [ ] 7.2 Interpolate global results onto the submodel boundary; driven BCs/loads; submodel solve

## 8. High-cycle fatigue (spec: high-cycle-fatigue)

- [ ] 8.1 `*HCF` evaluation over prior dynamic/modal results; worst-case location/life
- [ ] 8.2 Require a valid preceding results source (error if absent)

## 9. Parser, bindings, validation

- [ ] 9.1 Extend the deck parser to the Phase-5 card set (completes the recognized card universe)
- [ ] 9.2 Extend Python bindings + API-parity coverage test (full parity achieved)
- [ ] 9.3 Phase-5 reference corpus (CFD, network, EM, crack, optimization, refinement, submodel, HCF) + tolerances
- [ ] 9.4 pytest regression: each Phase-5 deck matches reference CalculiX within tolerance
- [ ] 9.5 Update README/docs; confirm full end-to-end coverage across all phases
