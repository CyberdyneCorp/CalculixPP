# Design — add-geometric-nonlinearity

## Context

`*BUCKLE` (backlog 2.2) is the last unshipped Phase-4 eigen procedure. Its only missing enabler is
the geometric (initial-stress) stiffness `K_geo`. This change scopes strictly to the **linear**
`K_geo` about a prestress state and wires `*BUCKLE` end-to-end as its consumer. Full finite-strain
NLGEOM is deferred: `*BUCKLE` needs only (1) a linear prestress solve for the reference
Gauss-point stress field, and (2) an initial-stress stiffness assembled **about the reference
configuration** — exactly CalculiX's `buckling=1` branch in `e_c3d.f`, which uses the block-diagonal
`senergy = ∇N_a·σ·∇N_b` term and **not** the material-tangent large-deformation branch.

## Approach

### Geometric (initial-stress) stiffness — reference-configuration linear buckling form

Reference: CalculiX `e_c3d.f` `buckling=1` branch. For element `e`, at each Gauss point `q` with
physical shape-function gradients `g_a = ∇N_a = dN_a/dX` (already computed by `physical_gradients`,
`element.cpp:361-380`, which returns `g` WITHOUT the Jacobian and scales by `det·w` separately) and
the reference Cauchy stress `σ` in Voigt6 `{σxx,σyy,σzz,σxy,σxz,σyz}` recovered from the prestress
solve, define the scalar

```
senergy(a,b) = g_a · σ · g_b
             = σxx·gax·gbx + σyy·gay·gby + σzz·gaz·gbz
               + σxy·(gax·gby + gay·gbx)
               + σxz·(gax·gbz + gaz·gbx)
               + σyz·(gay·gbz + gaz·gby)
```

The 3×3 nodal block of `K_geo` is **block-diagonal** in the spatial index `i` — the initial-stress
term couples the same translational direction of node `a` with node `b`:

```
K_geo[3a+i][3b+j] = δ_ij · Σ_gp senergy(a,b) · det(J) · w        (i,j ∈ {0,1,2})
```

This is the standard `Gᵀ Σ G` initial-stress matrix specialized to a symmetric stress and
small-strain (the off-diagonal-in-`i` terms vanish because `Σ = diag(σ,σ,σ)` in the 9×9 stress-block
sense). Element `K_geo` is `3n×3n`, symmetric; assemble exactly like `element_stiffness`
(`Σ_gp ... det·w`, upper triangle + `mirror_upper`).

**CalculiX-fidelity check (verified):** CalculiX `w(i,j) = shpj(i,a)·shpj(j,b)` with
`shpj = dN/dx·√xsj`, so `senergy = (s11·w11 + s12·(w12+w21) + … + s33·w33)·weight =
∇N_a·σ·∇N_b·detJ·w`, placed block-diagonal on the 3 translational DOFs. CalculiX puts `K` into `s`
(positive) and `K_geo` into `sm` with a **negative** sign (`sm = sm − senergyb`), giving the pencil
`K φ = (1/λ)(−K_geo) φ` — the orientation used below. No Jacobian double-count: CalculiX folds
`√xsj` into both `shpj`, while CalculixPP's `physical_gradients` returns `g` and scales by `det·w`
separately — algebraically identical.

### Two-step prestress driver

- **Step A (prestress):** solve the linear static problem `K u0 = f_ref` for the step's reference
  load `f_ref` (the `*BUCKLE` step's `*CLOAD`/`*DLOAD`). Recover Gauss-point stress `σ_q^e` via the
  shipped stress path. This IS the existing `solve_linear_static` + recover machinery.
- **Step B (eigenproblem):** assemble `K_geo(σ_ref)` and solve the generalized symmetric
  eigenproblem

  ```
  (K + λ K_geo) φ = 0    ⟺    K φ = μ (−K_geo) φ,   μ = 1/λ
  ```

  Solved as `(−K_geo) φ = θ K φ` with `θ = −1/λ` and `K` SPD (positive-definite for a constrained
  model). The critical buckling factor `λ_crit` is the smallest `λ > 0` — the load multiplier at
  which `K + λ K_geo` first becomes singular; report `λ_i` ascending-positive with mode shapes.
  Critical load = `λ_crit · f_ref`.

**SPD anchor.** The Cholesky anchor is the **unloaded** `K` (SPD before the reference load is
applied), which is exactly the classical linear-buckling formulation and keeps the reduction
`Â = L⁻¹(−K_geo)L⁻ᵀ` well-posed while `−K_geo` is indefinite.

## Eigen-extractor — a DEDICATED buckling routine (not literal reuse of `dense_extract_modes`)

The shipped `dense_extract_modes(K, M)` Choleskys its **second** argument and solves
`arg1 x = λ arg2 x`, then **clamps eigenvalues to ≥ 0** and fills `omega`/`frequency`. For buckling
that clamp and those fields are wrong. So `extract_buckling_modes(K, Kgeo, n)` is a **new**
extractor:

- Invoke the generalized reduction as `A = −K_geo`, `B = K` (SPD): Cholesky `K = LLᵀ`, form
  `Â = L⁻¹(−K_geo)L⁻ᵀ`, `eigh(Â) → θ`, back-transform `φ = L⁻ᵀ z`.
- Map `λ = −1/θ`; **filter to positive λ** (reject/skip non-positive or `θ→0` rigid-body / near-null
  `K_geo` directions, which map to `λ→±∞`); sort ascending; expand shapes through the constraint
  transform. No `≥0` clamp; the `eigenvalue` field holds the buckling factor λ (not `omega²`).

**Primary path = DENSE**, correct and validatable at `beamb` scale immediately. The scalable sparse
path is wired once the upstream SciPP target-selection ask lands (below).

## Stress-recovery refactor (regression-locked)

`strain_stress_at` lives in an **anonymous namespace** in `stress.cpp` (not exported).
`recover_gauss_stress(model, u) → vector<vector<Voigt6>>` (per-element per-Gauss reference stress)
must be factored out **without changing `recover_fields` byte output** — a real refactor with its
own regression-lock task, not a trivial add. Thermal correction is all-zero on the mechanical
buckling path. `recover_fields` delegates to / shares the extracted helper and stays byte-identical.

## Files to touch

- `include/calculixpp/fem/element.hpp` / `src/fem/element.cpp` — declare + implement
  `element_geometric_stiffness(type, coords, gp_stress)` returning the `3n×3n` block-diagonal
  `K_geo`; add an `accumulate_kgeo` helper next to `accumulate_ktangent`/`accumulate_fint`.
- `include/calculixpp/fem/assembly.hpp` / `src/fem/assembly.cpp` — declare + implement
  `assemble_geometric_stiffness(model, gp_stress)` mirroring `assemble_linear_static`
  (`build_dof_map`, gather coords, element kernel, `scatter_tangent` through the transform,
  `flush_coo`); respect `element_active_mask`. A zero stress field yields a zero matrix.
- `include/calculixpp/fem/stress.hpp` / `src/fem/stress.cpp` — add `recover_gauss_stress` (factor
  the anonymous-namespace `strain_stress_at`); keep `recover_fields` byte-identical.
- `include/calculixpp/numerics/eigensolution.hpp` / `src/numerics/eigensolution.cpp` — add the
  dedicated `extract_buckling_modes(K, Kgeo, n)` on the `(−K_geo, K)` pencil with `λ = −1/θ`.
- `include/calculixpp/numerics/buckling.hpp` / `src/numerics/buckling.cpp` (new) — `BucklingReport
  {factors; …}` and `solve_buckling(model, num_modes)` wiring the two steps, keeping the static
  drivers untouched.
- `include/calculixpp/core/model.hpp` — `Procedure::Buckling` (after `Frequency`) +
  `int num_buckling_modes{0}`.
- `src/io/inp_parser.cpp` — remove `*BUCKLE` from the deferred map; `begin_buckle()` sets
  `Procedure::Buckling`; `buckle_data()` reads the mode count (first field) into
  `num_buckling_modes`; accept and ignore the ARPACK-style tol/ncv/maxiter fields.
- `apps/cli/main.cpp` — a `Procedure::Buckling` dispatch branch (mirroring `Frequency`) calling
  `solve_buckling` and printing the `BUCKLING FACTOR OUTPUT` table (`MODE NO` / `BUCKLING FACTOR`).
- `python/bindings.cpp` — `buckling_result_dict` + a `Procedure::Buckling` case returning
  `{factors, mode_shapes, backend}`; add the name to the summary ladder.
- `src/io/writer/*.cpp` — `.dat` `BUCKLING FACTOR OUTPUT` block + `.frd` buckling mode datasets
  (**scoped as a thin follow-on** — see Scope).

## Decisions

1. **New `geometric-stiffness` capability, not folded into static-analysis.** `K_geo` is an operator
   reused by 3 distinct procedures (`*BUCKLE`, later `*STATIC, PERTURBATION`, large-strain),
   mirroring `element_mass`/`assemble_mass` under the eigen work. Static-analysis stays about the
   linear/Newton pipeline.
2. **Dedicated `assemble_geometric_stiffness` kernel, not the material-point path.** The material
   tangent path (`assemble_material_tangent`) is for finite-strain follow-on work. For `*BUCKLE` the
   lowest-risk route is a standalone kernel fed by recovered stresses, reusing the shipped, validated
   scatter/transform machinery and keeping the linear path byte-identical.
3. **Dedicated buckling extractor, not `dense_extract_modes` reuse** — the `≥0` clamp and
   `omega`/`frequency` fields are wrong for buckling factors (reviewer correction).

## Alternatives considered and rejected

1. **Full NLGEOM now** — balloons scope, not needed for `*BUCKLE`; the backlog explicitly allows it
   as a follow-on. Hold the line: keep the material interface and `strain_from_gradients` unchanged.
2. **Insert `K_geo` inside `element_tangent_force` behind a flag** — entangles the Newton path and
   needs per-Gauss stress plumbed through the material interface. A standalone, independently
   testable kernel is cleaner.
3. **Cast to `K x = μ(−K_geo)x` and call `eigsh`** — `eigsh` REQUIRES the M-argument SPD and
   `−K_geo` is indefinite. The correct pencil is `(A = −K_geo, B = K)` with `B = K` SPD, but
   `eigsh`'s `sigma`/`which` semantics return eigenvalues nearest `sigma` — the WRONG end for the
   critical (smallest positive) buckling factor. This is the gating upstream issue.

## Upstream asks

- **SciPP eigsh — target-selection for the buckling pencil.** `scipp::sparse::eigsh(K, M, k, sigma,
  …)` requires `M` SPD and (per its shift-invert Lanczos, factoring `(K−σM)` and returning
  eigenvalues nearest a shift) returns eigenvalues near `sigma`. For buckling the SPD matrix is `K`
  (not `K_geo`), and we need the SMALLEST POSITIVE `λ = −1/θ` of `(−K_geo)φ = θ K φ`, which is NOT
  "nearest `sigma=0`". Filed as **SciPP#18** (companion to SciPP#12/#15) — asks SciPP to either (a) add a `which`/mode selection
  targeting the smallest positive λ of a `K + λK_geo` buckling pencil, or (b) expose a generalized
  `eigsh` variant taking the SPD matrix as the `B`/mass argument with an indefinite `A` stiffness,
  returning the θ nearest a caller-supplied target so we can walk to the critical factor. Until this
  lands, the DENSE generalized path is used — correct and validatable but `O(n³)`; large models
  (>~1500 free DOF) have no scalable path.
- **(Follow-on) inertia / Sturm-sequence count** from the `K + λK_geo` factorization for missed-mode
  verification — tracked as backlog 1.3, becomes relevant once buckling ships; not required here.

## Fidelity gap (documented, acceptable)

CalculiX uses an ADAPTIVE `sigma` strategy (`arpackbu.c:504-508`: if the buckling factor is outside
`[5·sigma, 50000·sigma]`, reset `sigma = factor/500` and re-solve, via ARPACK buckling mode 4:
`iparam[6]=4`, `which=LM`, `sigma=1.0`). The dense/eigsh approach here does **not** replicate this
shift-adaptation — acceptable for the dense oracle, but it is precisely the target-selection
capability the upstream SciPP ask (**SciPP#18**) must supply for the sparse path.

## Scope note

This change sits at the upper edge for one OpenSpec change (new element kernel, new assembly kernel,
stress-recovery refactor, new eigen-extractor, new numerics driver, parser+enum, CLI, bindings,
writers). The coherent, archivable core is: `K_geo` kernel + assembly + prestress driver + dense
extractor + parser + one CLI/validation path. If the change grows too large to archive cleanly, the
CalculiX-exact `.dat`/`.frd` formatting (writer tasks) and half of the Python bindings are the first
to split into a thin output-plumbing follow-on — they are not the enabler math and CalculiX-exact
`.dat` formatting is a known time-sink. They are kept here but flagged.

## Risks

- **Eigensolver target-selection** — with `eigsh` unable to target the smallest positive factor, the
  sparse path is blocked; mitigated by shipping the dense generalized path (works at `beamb` scale)
  and filing the SciPP ask. Large models stay `O(n³)` until upstream lands.
- **Indefinite operator numerics** — `−K_geo` is indefinite and `K` can be near-singular under
  near-critical prestress. The reduction `Â = L⁻¹(−K_geo)L⁻ᵀ` (`K = LLᵀ`) is well-posed only while
  `K` is SPD; guard by anchoring the Cholesky on the **unloaded** `K` (SPD pre-load) — the classical
  linear-buckling formulation.
- **Sign/branch convention** — mapping `λ = −1/θ` and filtering to positive ascending factors is
  error-prone; validate the sign against `beamb.dat.ref` (all-positive factors) and reject spurious
  non-positive / infinite `θ` (rigid-body / near-null `K_geo` directions).
- **Element-family coverage** — `beamb.inp` is C3D20R (reduced integration). `K_geo` MUST use the
  SAME Gauss rule as the stiffness for consistency; verify reduced integration and C3D8R hourglass
  interaction do not corrupt `K_geo` (`K_geo` has no hourglass term — only the stress-gradient
  product).
- **Scope creep into finite-strain** — keep the material interface and `strain_from_gradients`
  UNCHANGED; hold `*STATIC, PERTURBATION` and large-strain as follow-on rows.
- **Prestress recovery reuse** — `recover_gauss_stress` must not alter `recover_fields` output;
  regression-lock the existing stress path.

## Validation vs CalculiX

- **Primary end-to-end** — run the stock `/home/leonardo/work/CalculiX/test/beamb.inp` (C3D20R Euler
  column, `*BUCKLE 10, 1.e-2`, unit reference load) through `*BUCKLE`; compare the lowest 10 factors
  to `beamb.dat.ref` (λ_1 = 0.4815456E+02 = 48.15, λ_2 = 0.1063175E+03 = 106.3, …) at a documented
  relative tolerance (target rel-L2 ~1e-3, consistent with the beam8f frequency validation). Also
  `beamb2.inp`.
- **Analytical Euler column** — mesh a slender prismatic bar (known `E`, `I`, `L`, BCs) with C3D8 and
  C3D20, apply a unit axial compressive reference load, and check `λ_crit·P_ref == π²EI/(kL)²` for
  pinned-pinned (`k=1`) and clamped-free (`k=2`) — a mesh-independent physics check.
- **Operator unit tests** — `element_geometric_stiffness` on one C3D8 under a prescribed uniaxial
  `σxx` equals the analytic block-diagonal `∇N·σ·∇N` matrix; `assemble_geometric_stiffness` shares
  `dof_eq`/`n_free` with `assemble_linear_static` and is symmetric; zero-stress → zero matrix.
- **Regression** — assert the linear-static, Newton, and `recover_fields` outputs are byte-identical
  after factoring `recover_gauss_stress` (no-op on non-buckling decks).
- **Python end-to-end** — solve `beamb.inp` via the bindings and assert the returned factors match
  the CLI/`.dat`, closing the loop through the Python API per the eigensolution reference-fidelity
  requirement.
