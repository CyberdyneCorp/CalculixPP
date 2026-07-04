# Tasks — add-submodeling

## 1. Model + core data

- [x] 1.1 Add `GlobalSolution` (global mesh view + per-node displacement) to a new
  `core/submodel.hpp`: node ids/coords, element type + connectivity, `displacement` aligned
  with the global nodes.
- [x] 1.2 Add `SubmodelSpec` (boundary node set name, global elset name, `TYPE=NODE`) and a
  `std::vector<SubmodelSpec> submodels` + a `driven` flag on `Spc` to the model.
- [x] 1.3 Add `Procedure::Submodel`? — NOT needed; a submodel is an ordinary `*STATIC` step whose
  driven SPCs are filled before the solve. Keep `Procedure::Static`; gate on `!submodels.empty()`.

## 2. Interpolation kernel (numerics/submodel.{hpp,cpp})

- [x] 2.1 `natural_coords(elem, coords, X)` — Newton inverse-isoparametric map: solve
  `x(ξ) = X` for the natural coords ξ of a physical point X inside a global element, using
  `fem::shape` for `x(ξ)` and `physical_gradients`/the Jacobian for `dx/dξ`. Return ξ + the
  residual distance.
- [x] 2.2 `find_host_element(global, X)` — search the global elements for the one containing X
  (ξ inside the reference domain within tolerance / smallest residual). Throw if none is found.
- [x] 2.3 `interpolate_displacement(global, X)` — locate the host element, evaluate the shape
  functions at ξ, and return `Σ N_a(ξ) U_a` (the interpolated global displacement at X).
- [x] 2.4 `solve_submodel(model, global)` — for every driven `*BOUNDARY, SUBMODEL` node, set its
  SPC value to the interpolated global displacement component, then run `solve_linear_static`.

## 3. Parser

- [x] 3.1 `*SUBMODEL` card: `TYPE=NODE` (+ global elset) records a `SubmodelSpec`; register the
  boundary node set as driven. Remove `*SUBMODEL` from the deferred/rejected cards.
- [x] 3.2 `*BOUNDARY, SUBMODEL`: mark the emitted SPCs as `driven` (value filled at solve time).
- [x] 3.3 Reject `*SUBMODEL, TYPE=SURFACE` and `*CLOAD/*DSLOAD, SUBMODEL` with a clear
  "deferred" message (never silently ignore). [~]

## 4. Bindings + CLI

- [x] 4.1 Python `solve_submodel(sub_text, global_ids, global_coords, global_conn, global_type,
  global_disp, solver, backend)` — build a `GlobalSolution`, parse the submodel deck, run
  `solve_submodel`, return the standard result dict (displacement/stress/reaction).
- [x] 4.2 CLI: `*SUBMODEL` decks need a global result file on disk (deferred), so the CLI reports a
  clear "provide the global solution via the Python API" message rather than mis-solving. [~]

## 5. Validation + regression

- [x] 5.1 `tests/test_submodel.cpp`: build a coarse global beam deck + solution, cut a sub-region,
  drive its cut-boundary nodes from the interpolated global displacements, solve, and assert the
  submodel reproduces the global displacement inside the region to rel-L2 < 1e-3. Also assert a
  node outside the global element set throws.
- [x] 5.2 Wire `test_submodel` into `tests/CMakeLists.txt`; run `ctest` and the Python tests.

## Deferred (with reasons)

- **`.frd`/results-file global source** [~] — the interpolation kernel is source-agnostic; reading
  and indexing a global `.frd` from disk is I/O plumbing orthogonal to the algorithm. The in-memory
  `GlobalSolution` (built in tests/bindings) exercises the full interpolation + solve path.
- **`TYPE=SURFACE` boundary** [~] — the same host-element interpolation applies to surface nodes;
  only the boundary-entity enumeration differs. Deferred to keep the slice to node sets.
- **Force/pressure-driven `*CLOAD/*DSLOAD, SUBMODEL` and temperature-driven submodels** [~] —
  need interpolation of the global stress/temperature field (extrapolated to nodes) rather than the
  primary displacement field; a distinct follow-on. The displacement-driven case is by far the most
  common and is the one validated here.
