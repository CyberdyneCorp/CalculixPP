# Tasks — add-geometric-nonlinearity

Implements the linear geometric (initial-stress) stiffness enabler and wires `*BUCKLE`
end-to-end as its consumer. Full finite-strain NLGEOM is out of scope (follow-on rows). Depends on
`phase-1-foundation` (assembly/solve spine, stress recovery) and
`phase-4-dynamics-and-eigenproblems` (eigensolution engine).

## 1. Geometric-stiffness element kernel (spec: geometric-stiffness [NEW])

- [ ] 1.1 Declare `element_geometric_stiffness(ElementType, std::span<const Vec3> coords, const
      std::vector<Voigt6>& gp_stress)` in `include/calculixpp/fem/element.hpp`, returning the
      `3n×3n` block-diagonal `K_geo` (one `Voigt6` reference stress per Gauss point); document the
      `senergy` formula and reference-config initial-stress scope, mirroring `element_mass` /
      `element_stiffness`.
- [ ] 1.2 Implement `element_geometric_stiffness` in `src/fem/element.cpp`: loop the SAME
      `gauss_rule` as the stiffness, `physical_gradients → g`, accumulate `senergy(a,b)·det·w` into
      the 3 diagonal (`i==j`) DOF pairs (upper triangle + `mirror_upper`); add an `accumulate_kgeo`
      helper next to `accumulate_ktangent`/`accumulate_fint`.
- [ ] 1.3 Unit-test on a single C3D8 under a prescribed uniaxial `σxx` that the block-diagonal
      `senergy` assembly reproduces the analytic `∇N·σ·∇N` initial-stress matrix, and that a zero
      stress field gives a zero element matrix.

## 2. Stress-recovery refactor (spec: geometric-stiffness [NEW])

- [ ] 2.1 Add `recover_gauss_stress(const Model&, const std::vector<Vec3>& u) →
      std::vector<std::vector<Voigt6>>` (per-element per-Gauss reference stress) in
      `include/calculixpp/fem/stress.hpp`, factoring the anonymous-namespace `strain_stress_at` out
      of `src/fem/stress.cpp` (thermal correction all-zero on the mechanical buckling path).
- [ ] 2.2 **[regression-lock]** Assert `recover_fields` output is byte-identical after the
      extraction (`recover_fields` delegates to / shares the helper). Blocks 5.x.

## 3. Geometric-stiffness assembly kernel (spec: geometric-stiffness [NEW])

- [ ] 3.1 Declare `assemble_geometric_stiffness(const Model&, const
      std::vector<std::vector<Voigt6>>& gp_stress) → LinearSystem` in
      `include/calculixpp/fem/assembly.hpp`; document same free-DOF numbering / constraint transform
      as `assemble_linear_static`/`assemble_mass`, and that a zero stress field yields a zero matrix.
- [ ] 3.2 Implement it in `src/fem/assembly.cpp` mirroring `assemble_linear_static` (`build_dof_map`,
      gather coords, `element_geometric_stiffness`, `scatter_tangent` through the transform,
      `flush_coo`); respect `element_active_mask`.
- [ ] 3.3 Test that its COO free-DOF count and `dof_eq` match `assemble_linear_static` for the same
      model, that it is symmetric, and that a zero stress field yields an all-zero matrix.

## 4. Buckling eigen-extractor (spec: modal-and-buckling-analysis [MODIFIED])

- [ ] 4.1 **[upstream-gated for the SPARSE path]** File the SciPP target-selection ask (companion to
      SciPP#15): buckling pencil `(A = −K_geo, B = K)`, need smallest positive `λ = −1/θ`.
- [ ] 4.2 Add a DEDICATED `extract_buckling_modes(const fem::LinearSystem& K, const fem::LinearSystem&
      Kgeo, std::size_t num_modes) → EigenBasis` in `include/calculixpp/numerics/eigensolution.hpp`
      whose `eigenvalue` field holds the buckling factor λ (ascending positive) and `shape` the
      buckling mode; document the `(−K_geo, K)` pencil and `λ = −1/θ` mapping. **Not** a reuse of
      `dense_extract_modes` — the `≥0` clamp and `omega`/`frequency` fields are wrong for buckling.
- [ ] 4.3 Implement it in `src/numerics/eigensolution.cpp` via the DENSE generalized path first
      (Cholesky of the SPD **unloaded** `K` + `eigh` on `Â = L⁻¹(−K_geo)L⁻ᵀ`): map `λ = −1/θ`,
      filter to positive λ (reject non-positive / `θ→0` rigid-body directions), sort ascending,
      normalize, expand shapes through the constraint transform. Wire the sparse `eigsh` path once
      4.1 lands.

## 5. Two-step prestress driver (spec: geometric-stiffness [NEW])

- [ ] 5.1 Add `struct BucklingReport { std::vector<Real> factors; … }` and
      `solve_buckling(const Model&, std::size_t num_modes)` in a NEW
      `include/calculixpp/numerics/buckling.hpp` / `src/numerics/buckling.cpp` (keep static drivers
      untouched): Step A `solve_linear_static` prestress → `recover_gauss_stress` → Step B assemble
      `K` and `K_geo` → `extract_buckling_modes` → ascending-positive factors + shapes. Requires 2.2,
      3.x, 4.x.
- [ ] 5.2 Assert the linear-static and Newton paths are byte-for-byte unchanged for a deck with no
      `*BUCKLE`/geometric step.

## 6. Parser + model (spec: input-deck-parsing [MODIFIED])

- [ ] 6.1 Add `Procedure::Buckling` (after `Frequency`) and `int num_buckling_modes{0}` to
      `include/calculixpp/core/model.hpp`.
- [ ] 6.2 In `src/io/inp_parser.cpp`: remove `*BUCKLE` from the deferred-card map; add
      `begin_buckle()` (sets `Procedure::Buckling`), `buckle_data()` (reads the requested mode count
      from the first field into `num_buckling_modes`, accepts + ignores the ARPACK-style
      tol/ncv/maxiter fields), and register both in the card + data dispatch tables.
- [ ] 6.3 Parser test: `beamb.inp` parses to `Procedure::Buckling` with `num_buckling_modes == 10`;
      the trailing `1.e-2` tolerance field is accepted and ignored without error.

## 7. CLI + Python wiring (spec: modal-and-buckling-analysis [MODIFIED])

- [ ] 7.1 Add a `Procedure::Buckling` dispatch branch in `apps/cli/main.cpp` (mirroring `Frequency`)
      calling `solve_buckling` and printing a `BUCKLING FACTOR OUTPUT` table (`MODE NO` /
      `BUCKLING FACTOR`).
- [ ] 7.2 Add `buckling_result_dict` and a `Procedure::Buckling` case in `python/bindings.cpp`
      returning `{factors, mode_shapes, backend}`; add the procedure name to the summary string
      ladder. **[may split to a thin follow-on if the change grows too large to archive cleanly.]**

## 8. Result writers (spec: modal-and-buckling-analysis [MODIFIED])

- [ ] 8.1 Emit the `.dat` `BUCKLING FACTOR OUTPUT` block and the `.frd` buckling mode datasets in
      `src/io/writer/*.cpp`, matching `beamb.dat.ref` formatting for the reference comparison.
      **[output plumbing; CalculiX-exact `.dat` formatting is a known time-sink — first to split to
      a thin follow-on if archiving pressure demands it.]**

## 9. Validation vs CalculiX (spec: modal-and-buckling-analysis [MODIFIED])

- [ ] 9.1 Run `beamb.inp` end-to-end; compare λ_1..λ_10 to `beamb.dat.ref`
      (λ_1 = 0.4815456E+02 = 48.15, λ_2 = 0.1063175E+03 = 106.3, …) within a documented relative
      tolerance (target rel-L2 ~1e-3). Also run `beamb2.inp`.
- [ ] 9.2 Analytical Euler-column regression: mesh a slender prismatic bar with C3D8 and C3D20 under
      a unit axial compressive reference load; check `λ_crit·P_ref == π²EI/(kL)²` for pinned-pinned
      (`k=1`) and clamped-free (`k=2`) — mesh-independent physics.
- [ ] 9.3 Python end-to-end: solve `beamb.inp` via the bindings and assert the returned factors match
      the CLI/`.dat` output.

## 10. OpenSpec + backlog hygiene

- [ ] 10.1 On archive, mark the linear-`K_geo` slice of backlog row **2.2 `*BUCKLE`** done; leave
      finite-strain NLGEOM + `*STATIC, PERTURBATION` + large-strain as follow-on rows in `BACKLOG.md`
      (with the NLGEOM cross-cutting-enabler note pointing at this change for the linear slice).
