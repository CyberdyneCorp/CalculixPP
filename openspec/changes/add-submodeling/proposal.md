## Why

The `submodeling` capability is specified but unimplemented: the parser rejects `*SUBMODEL`
and there is no driver that reads a global solution, interpolates it onto a finer local mesh,
and solves the local model. Submodeling is the standard workflow for resolving a stress
concentration in a small region at high mesh density without re-solving the whole part — the
local model is driven purely by the global displacement field on its cut boundary.

This change delivers the smallest correct end-to-end slice: the **displacement-driven** case.
A `*SUBMODEL, TYPE=NODE` card names the driven boundary node set and the global result source;
`*BOUNDARY, SUBMODEL` prescribes each listed DOF of those nodes to the global displacement
interpolated at the node's location; the submodel then runs as an ordinary `*STATIC` step.

## What Changes

- **Global-result interpolation kernel** — a `numerics/submodel.{hpp,cpp}` driver that, for each
  submodel boundary node, locates the host global element containing the node's coordinates,
  computes the node's natural coordinates by a Newton inverse-isoparametric map, and interpolates
  the global nodal displacement field with that element's shape functions. It reuses the existing
  `fem::shape` / `physical_gradients` machinery — no new element math. A node outside every host
  element throws a clear error.
- **`GlobalSolution`** — a lightweight in-memory carrier for the global mesh + nodal displacement
  field (node coords, element connectivity/type, per-node U). Built directly in tests/bindings;
  reading a global `.frd` from disk is a thin follow-on (see deferred).
- **`*SUBMODEL` parser card** — `TYPE=NODE` records the driven boundary node set and the global
  element set searched for host elements. Removed from the deferred-card map.
- **`*BOUNDARY, SUBMODEL`** — the `SUBMODEL` parameter flags a boundary card as driven: its DOFs
  are prescribed to interpolated global displacements at solve time rather than to literal values.
- **Submodel solve driver + dispatch** — `solve_submodel(model, global)` interpolates the driven
  values into the model's SPCs and runs the linear-static solve. Reachable from the Python bindings
  (`solve_submodel`) returning the standard displacement/stress/reaction result dict.

## Capabilities

### Modified Capabilities

- `submodeling`: the baseline requirements are refined to the implemented displacement-driven slice
  — `TYPE=NODE` boundary declaration, host-element location + shape-function interpolation of the
  global displacement field, `*BOUNDARY, SUBMODEL` prescribing the interpolated displacements, the
  ordinary local `*STATIC` solve through the compute backend, and the cut-from-a-solved-beam
  fidelity oracle. Stress/temperature-driven submodels (`*CLOAD/*DSLOAD, SUBMODEL`, `TYPE=SURFACE`)
  are scoped as deferred with a note.

## Impact

- **New code:** `numerics/submodel.{hpp,cpp}` (host-element search, inverse-isoparametric Newton
  map, displacement interpolation, `solve_submodel` driver); `GlobalSolution` + `SubmodelSpec` on
  the model; `*SUBMODEL` parser dispatch and the `SUBMODEL` flag on `*BOUNDARY`.
- **Extended code:** Python binding `solve_submodel(...)` returning the standard result dict; a
  `test_submodel` C++ unit test building the global + submodel decks inline.
- **Validation oracle:** a beam solved globally (coarse), then a sub-region cut out and re-solved
  driven by the interpolated global displacements on its cut faces — the submodel displacement
  field reproduces the global field inside the region to rel-L2 < 1e-3.
- **Dependencies:** unchanged (NumPP, SciPP, pybind11; no GPU required).
- **Deferred (out of scope, documented in tasks):** reading the global field from a `.frd`/results
  file on disk; `TYPE=SURFACE` boundaries; force/pressure-driven `*CLOAD/*DSLOAD, SUBMODEL`;
  temperature-driven submodels; cross-step global-step selection.
- **Depends on:** `phase-1-foundation` (mesh, shape functions, linear-static solve, constraint
  transform).
