# Tasks — add design-sensitivity core

## 1. Model + parser

- [x] 1.1 Add `Procedure::Sensitivity` to `core/model.hpp`.
- [x] 1.2 Add design data to `Model`: coordinate design variables (node id + component
      1..3) and a list of `DesignResponse` (name, kind = Compliance / StrainEnergy /
      Displacement, optional target node/dof for Displacement).
- [x] 1.3 Parse `*DESIGN VARIABLES, TYPE=COORDINATE` (data line: nset/node) — expand to
      per-node x/y/z coordinate variables.
- [x] 1.4 Parse `*DESIGN RESPONSE` / `*OBJECTIVE` (name + response type; `ALL-DISP` /
      nodal `U` displacement target). Accept `STRAIN ENERGY` (compliance) and a nodal
      displacement response.
- [x] 1.5 Parse the `*SENSITIVITY` step card → `Procedure::Sensitivity`.
- [~] 1.6 `*ORIENTATION` design variables (`sensi_orien`) — DEFERRED (only
      `TYPE=COORDINATE` shape variables this slice); card rejected with a clear message.

## 2. Sensitivity driver (`numerics/sensitivity`)

- [x] 2.1 Solve the primal `K u = f` (reuse `assemble_linear_static` + backend solve).
- [x] 2.2 Semi-analytic operator derivatives: for each coordinate design variable,
      perturb the coordinate by ±h, re-assemble ONLY the affected elements' stiffness /
      load, form dK/dx and df/dx by central difference of the element operators.
- [x] 2.3 Compliance (self-adjoint) sensitivity: `dg/dx = 2 uᵀ df/dx − uᵀ dK/dx u`
      (reuses `u`; no adjoint solve).
- [x] 2.4 Nodal-displacement (general linear) sensitivity: adjoint `K λ = c` (reuses the
      primal factorization / same operator), `dg/dx = λᵀ (df/dx − dK/dx u)`.
- [~] 2.5 `*FILTER` sensitivity smoothing — DEFERRED (regularization layer on top of the
      raw gradient; not needed for correctness).

## 3. Dispatch + bindings + test

- [x] 3.1 CLI: dispatch `Procedure::Sensitivity`, print the gradient table.
- [x] 3.2 Python `solve()` result dict: `procedure`, per-variable `gradient`,
      `design_var_node` / `design_var_comp`, `response_name`.
- [x] 3.3 C++ unit test `test_sensitivity`: **finite-difference gradient check** —
      adjoint dObjective/dx vs central `(O(+h)−O(−h))/2h` re-solve, <1e-4 relative, on a
      small cantilever beam, for BOTH compliance and a nodal-displacement response.

## 4. Deferred sub-features (documented, not faked)

- [~] 4.1 `*FEASIBLE DIRECTION` optimization loop — DEFERRED (consumes the gradient; a
      separate optimizer slice).
- [~] 4.2 Robust design `*ROBUST DESIGN` / `*RANDOM FIELD` / `*CORRELATION LENGTH` —
      DEFERRED (uncertainty layer; needs a random-field decomposition).

## 5. Validate + guard

- [x] 5.1 `openspec validate add-design-optimization --strict` passes.
- [x] 5.2 Clean build; `ctest` + python tests green (no regression).
