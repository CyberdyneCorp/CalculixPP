## Why

Phases 1–4 cover structural statics, thermal/contact, and dynamics — the structural core.
Phase 5 completes the solver: the remaining physics (3-D CFD and 1-D fluid networks,
electromagnetics, crack propagation, design sensitivity/optimization) and the
adaptivity/reuse capabilities (adaptive mesh refinement, submodeling, high-cycle fatigue).
The physics is already in the baseline; the one genuinely-new piece of machinery is the
multi-physics field-coupling engine that CFD, electromagnetic-thermal (Joule), and
network-thermal coupling all require.

## What Changes

- **CFD & 1-D networks** — 3-D compressible/incompressible Navier-Stokes (`*CFD`), 1-D
  fluid-network elements (orifices, pipes, valves, labyrinths), fluid/physical constants,
  face/flow boundary conditions, network MPCs, and coupling to heat transfer.
- **Electromagnetics** — magnetostatics, eddy-current electromagnetics, electric
  conduction (`*ELECTROMAGNETICS`), EM material properties, Joule-heating coupling to the
  thermal field, and the enveloping air domain.
- **Crack propagation** — stress-intensity factors along the crack front, crack length/shape
  advance with CyberCadKernel remeshing, growth-rate law, and crack output
  (`*CRACK PROPAGATION`).
- **Design sensitivity & optimization** — design variables (coordinate/orientation),
  responses/objectives/constraints, adjoint sensitivities, filtering, feasible-direction
  optimization, and robust design / random field.
- **Adaptive mesh refinement** — `*REFINE MESH` solve–refine–resolve loop with
  CyberCadKernel tet remeshing; `*INITIAL MESH` / `*USE REFINED MESH`.
- **Submodeling** — `*SUBMODEL` driven from a global run's boundary results.
- **High-cycle fatigue** — `*HCF` life/critical-location evaluation over prior dynamic/modal
  results.
- **New capability** — `field-coupling`: the segregated multi-physics coupling engine.
- No breaking changes. GPU backends remain optional. This is the final phase.

## Capabilities

### New Capabilities

- `field-coupling`: the reusable segregated multi-physics engine — orchestrate a coupled
  solve across fields (fluid momentum/pressure, temperature, electromagnetic potential),
  exchange coupling terms (Joule heat, buoyancy, forced-convection film, fluid tractions),
  and converge the coupled system. Peer to `linear-algebra-and-solvers`,
  `nonlinear-solution-control`, `contact-search`, and `eigensolution` (engine, not
  procedure); consumed by CFD, electromagnetic-thermal, and network-thermal coupling.

### Modified Capabilities

None. `cfd-and-network-analysis`, `electromagnetic-analysis`, `crack-propagation`,
`design-optimization`, `mesh-refinement`, `submodeling`, and `high-cycle-fatigue` already
exist as baseline capabilities; this change implements them. `tasks.md` cites the exact
baseline requirements. Only `field-coupling` is added as a delta. Adjoint sensitivity is
implemented within `design-optimization` and tet remeshing within `mesh-refinement` /
`crack-propagation` rather than extracted as separate engines this phase, to keep scope
bounded.

## Impact

- **New code:** the field-coupling engine (segregated iteration, coupling-term transfer,
  pressure-velocity coupling), CFD/network solvers, EM potential solvers + Joule coupling,
  crack-front SIF + remeshing, adjoint sensitivity + optimization drivers, adaptive
  refinement loop, submodel interpolation, and HCF post-processing.
- **Extended code:** parser coverage for the Phase-5 card set; Python bindings +
  API-parity check; the regression corpus grows with CFD, network, EM, crack, optimization,
  refinement, submodel, and HCF decks.
- **Dependencies:** unchanged (NumPP, CyberCadKernel, pybind11; no GPU required).
- **Depends on:** `phase-1-foundation` … `phase-4-dynamics-and-eigenproblems` (assembly/solve
  spine, nonlinear driver, thermal fields/contact search, eigensolution for
  frequency-based sensitivities and HCF).
