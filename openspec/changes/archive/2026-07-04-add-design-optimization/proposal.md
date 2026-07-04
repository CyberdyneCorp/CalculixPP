# Add design-sensitivity core (`*DESIGN VARIABLES` / `*DESIGN RESPONSE` / `*OBJECTIVE` / `*SENSITIVITY`)

## Why

The `design-optimization` capability is entirely unimplemented — every `*DESIGN
VARIABLES` / `*SENSITIVITY` deck currently fails with a deferral parse error. The
capability's hard, validatable core is the **design sensitivity** dObjective/dx of a
response with respect to coordinate design variables. That gradient is what every
gradient-based optimizer, filter, and feasible-direction step consumes, so it is the
right first slice: it is small, self-contained, and has an authoritative correctness
gate (a central finite-difference gradient check).

The adjoint method makes this cheap and exact: for the linear system `K(x) u = f(x)`,
the sensitivity of a linear/self-adjoint response **reuses the primal solution** (and,
for a general linear response, one extra solve with the SAME operator) rather than
re-solving per design variable. This is the CalculiX `sensi_coor` path
(`src/sensitivitys.f`, `objective_shapeener_tot.f`), reimplemented on our
assembly/ComputeBackend stack — never copied.

## What

- Parse `*DESIGN VARIABLES, TYPE=COORDINATE`, `*DESIGN RESPONSE`, `*OBJECTIVE`, and a
  `*SENSITIVITY` step; add `Procedure::Sensitivity` and the model data (design-variable
  node/component list, design responses).
- New `numerics/sensitivity` driver: solve the primal `K u = f`, then compute
  dObjective/dx for each coordinate design variable by the adjoint method, obtaining
  the operator derivatives dK/dx and df/dx by a **semi-analytic** perturbation of the
  element stiffness/load (the CalculiX approach), reusing the primal factorization for
  the adjoint solve.
- Supported responses this slice: **compliance** `fᵀu` (a.k.a. STRAIN ENERGY, the
  self-adjoint case) and a **single nodal displacement** DOF `cᵀu` (the general linear
  case). Both are exact linear-elastic responses with a closed-form adjoint.
- CLI dispatch, Python `solve()` result dict (`gradient` per design variable), and a
  C++ unit test whose gate is the finite-difference gradient check (<1e-4 relative).

### Deferred (parsed/marked, not faked)

- `*FILTER` sensitivity smoothing, the `*FEASIBLE DIRECTION` optimization loop, and
  robust design (`*ROBUST DESIGN` / `*RANDOM FIELD` / `*CORRELATION LENGTH`) — each is a
  layer *on top of* the sensitivity core and out of this slice. `*ORIENTATION` design
  variables (`sensi_orien`) are also deferred; only `TYPE=COORDINATE` is implemented.
- Mass / stress / eigenfrequency responses beyond compliance + nodal displacement.

## Impact

- Specs: `design-optimization` (refines "Design variables and responses" +
  "Sensitivity analysis" to the implemented behaviour).
- Code: `core/model.hpp` (Procedure + design data), `io/inp_parser.cpp`,
  `numerics/sensitivity.{hpp,cpp}`, `apps/cli/main.cpp`, `python/bindings.cpp`, a new
  C++ test. No change to existing analysis paths — a deck without `*DESIGN VARIABLES`
  is byte-for-byte unchanged.
