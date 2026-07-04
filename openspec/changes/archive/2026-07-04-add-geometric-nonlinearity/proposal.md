## Why

CalculiX++ ships the eigensolution engine, the linear-static solve, and the stress-recovery
machinery, but `*BUCKLE` (backlog item 2.2) is still unimplemented and the parser rejects the
card. The single missing enabler is a **geometric (initial-stress) stiffness operator** `K_geo`:
the reference-configuration matrix that, combined with the linear stiffness `K`, produces the
buckling pencil `(K + λ K_geo) φ = 0`. This is exactly the `buckling=1` branch of CalculiX
`e_c3d.f` — a linear prestress solve followed by an eigenproblem — and it does **not** require
full finite-strain kinematics.

This change delivers the smallest correct slice that unblocks the backlog: the linear `K_geo`
operator plus the two-step prestress driver, wired end-to-end so `*BUCKLE` runs and matches the
stock `beamb` reference (λ_1 = 48.15). Full finite-strain NLGEOM — updated deformation gradient
`F`, Green-Lagrange strain, PK2 stress, push-forward — is **explicitly deferred** to follow-on
rows (`*STATIC, PERTURBATION` and large-strain hyperelasticity), which `*BUCKLE` does not need.

## What Changes

- **New geometric-stiffness enabler** — an `assemble_geometric_stiffness(model, gp_stress)` kernel
  that builds the block-diagonal initial-stress matrix `K_geo[3a+i][3b+j] = δ_ij · Σ_gp (∇N_a · σ ·
  ∇N_b) · det(J) · w` from a supplied per-integration-point reference stress field, on the SAME
  free-DOF numbering and constraint transform as the linear stiffness. The linear-static and Newton
  paths are unchanged (a zero stress field yields a zero matrix).
- **Two-step prestress driver** — solve the linear static response to the step's reference load,
  recover the reference Gauss-point stresses, assemble `K` and `K_geo`, and solve the generalized
  eigenproblem. The stress recovery is factored out of the existing path with a byte-identical
  regression lock on `recover_fields`.
- **`*BUCKLE` made real** — the `modal-and-buckling-analysis` baseline `*BUCKLE` requirement is
  implemented and its scenarios refined; the buckling factors are returned ascending-positive with
  their mode shapes, reachable from the CLI and Python bindings.
- **Parser accepts `*BUCKLE`** — removed from the deferred-card map; the requested mode count is
  read from the data line and the ARPACK-style tolerance/subspace/iteration fields are accepted and
  ignored.
- **Buckling eigen-extractor** — a dedicated `extract_buckling_modes(K, Kgeo, n)` on the pencil
  `(A = −K_geo, B = K)` mapping `λ = −1/θ`, filtering to positive factors, sorted ascending. The
  DENSE generalized path (`numpp` Cholesky of the SPD `K` + `eigh`) is the primary, validatable
  path; the scalable sparse path is gated on an upstream SciPP target-selection ask (see design).

## Capabilities

### New Capabilities

- `geometric-stiffness`: the reusable `K_geo` operator (a new assembly kernel) plus the two-step
  prestress driver. It is an engine, not a procedure — cross-cutting and later consumed by
  `*STATIC, PERTURBATION` and large-strain follow-on work as well as by `*BUCKLE` — so it earns its
  own capability, mirroring how `element_mass`/`assemble_mass` sit under the eigen work rather than
  inside static-analysis. Keeping it separate keeps `static-analysis` about the linear/Newton
  pipeline.

### Modified Capabilities

- `modal-and-buckling-analysis`: the existing (unimplemented) `*BUCKLE` requirement is made real —
  refined to reference the `K_geo` prestress state, ascending-positive factors + mode shapes, CLI +
  Python reachability, and the `beamb` / analytical Euler-column validations.
- `input-deck-parsing`: a `*BUCKLE` procedure-card requirement is added — the parser accepts the
  card, records the buckling procedure and requested mode count, and ignores the trailing
  ARPACK-style fields, instead of rejecting it as deferred.

## Impact

- **New code:** the `element_geometric_stiffness` kernel and `assemble_geometric_stiffness`
  assembler, `recover_gauss_stress` (factored from the existing stress path), the
  `extract_buckling_modes` eigen-extractor, and a small `numerics/buckling.{hpp,cpp}` two-step
  driver. `Procedure::Buckling` + `num_buckling_modes` on the model; `*BUCKLE` parser dispatch.
- **Extended code:** CLI dispatch printing a `BUCKLING FACTOR OUTPUT` table and Python bindings
  returning the factors + mode shapes (output plumbing; the CalculiX-exact `.dat`/`.frd` formatting
  is scoped as a thin follow-on — see design, Scope).
- **Regression corpus:** the stock `beamb.inp` / `beamb2.inp` buckling decks and an analytical
  Euler-column check (pinned-pinned k=1, clamped-free k=2) meshed with C3D8/C3D20.
- **Dependencies:** unchanged (NumPP, SciPP, pybind11; no GPU required). One upstream SciPP
  target-selection ask filed to make the sparse buckling path scalable; the dense path ships now.
- **Deferred (out of scope):** finite-strain NLGEOM (`F`, Green-Lagrange, PK2, push-forward),
  `*STATIC, PERTURBATION`, large-strain hyperelasticity, and the CalculiX adaptive-`sigma`
  buckling-shift strategy — all follow-on rows.
- **Depends on:** `phase-1-foundation` (assembly/solve spine, stress recovery) and
  `phase-4-dynamics-and-eigenproblems` (the eigensolution engine).
