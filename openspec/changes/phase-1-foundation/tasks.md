# Phase 1 — Foundation & Linear-Static Vertical Slice

Tasks reference the capability they implement. Spec deltas exist only for
`build-and-tooling`; all other references are to existing baseline specs.

## 1. Repository & build scaffolding (spec: build-and-tooling)

- [ ] 1.1 Create module layout: `core/`, `numerics/`, `geometry/`, `cli/`, `python/`, `tests/`, `cmake/`
- [ ] 1.2 Top-level `CMakeLists.txt` enforcing C++20; core library + CLI targets
- [x] 1.3 Add `.clang-format`/`.clang-tidy` and wire warnings-as-errors on the core
- [ ] 1.4 README quickstart (configure, build, run a deck, run Python)

## 2. Dependency integration (spec: build-and-tooling)

- [x] 2.1 Integrate NumPP (`find_package`) + SciPP (`add_subdirectory`; provides the sparse module); smoke test a sparse solve. NOTE: sparse lives in SciPP, not NumPP; see `scripts/bootstrap_deps.sh`
- [ ] 2.2 Integrate CyberCadKernel at a pinned version; smoke test loading a B-rep
- [ ] 2.3 Integrate pybind11; build an empty `calculixpp` Python module
- [ ] 2.4 Clear configure-time failure when a dependency/version is unresolved

## 3. Compute backend — CPU (spec: compute-backend, linear-algebra-and-solvers)

- [ ] 3.1 Define `ComputeBackend` interface: assembly scatter, SpMV, sparse factor/solve
- [ ] 3.2 Implement CPU reference backend over SciPP sparse (+ NumPP dense/ndarray) (default, always available)
- [ ] 3.3 Backend selection surface (explicit → best available → CPU fallback), CPU-only for now
- [ ] 3.4 Guarantee build + tests pass with no GPU toolkit present (CI job)

## 4. Core FE data model (spec: mesh-and-model)

- [x] 4.1 Node store (ids, coordinates) and element store (connectivity, type)
- [x] 4.2 `C3D4` and `C3D10` tetrahedral element types
- [ ] 4.3 Node sets / element sets (`*NSET`/`*ELSET`) and surfaces (`*SURFACE`)
- [x] 4.4 DOF numbering / equation map
- [x] 4.5 Central `cpp_index_t` (32-bit signed) typedef used by all node/element/DOF indexing

## 5. Mesh import bypass (spec: mesh-processing)

- [ ] 5.1 Import mesh directly from the deck (`*NODE`/`*ELEMENT`), bypassing CAD meshing
- [ ] 5.2 Stub CyberCadKernel geometry-import/triangulate/tet interfaces (validated later phase)

## 6. Input deck parser — Phase-1 subset (spec: input-deck-parsing)

- [ ] 6.1 Keyword tokenizer/lexer for Abaqus-style cards (case-insensitive, comments, continuations)
- [ ] 6.2 Parse `*NODE`, `*ELEMENT`, `*NSET`, `*ELSET`, `*SURFACE`
- [ ] 6.3 Parse `*MATERIAL`, `*ELASTIC`, `*DENSITY`, `*SOLID SECTION`
- [ ] 6.4 Parse `*BOUNDARY`, `*CLOAD`, `*DLOAD` (pressure), `*STEP`/`*STATIC`, `*AMPLITUDE` (constant/step)
- [ ] 6.5 Actionable parse errors with card + line context

## 7. Material & section (spec: material-models, element-sections)

- [x] 7.1 Isotropic linear-elastic constitutive model (E, ν) and density
- [x] 7.2 `*SOLID SECTION` binding material to element set; local orientation hook (identity for now)

## 8. Assembly & loads (spec: static-analysis, loads-and-boundary-conditions)

- [x] 8.1 Element stiffness kernels for `C3D4`/`C3D10` (Gauss integration)
- [x] 8.2 Global assembly into COO triplets → SciPP `CsrMatrix::from_coo` (compute-backend abstraction still to be formalized)
- [~] 8.3 Concentrated loads (`*CLOAD`) done; pressure loads (`*DLOAD`) pending
- [x] 8.4 Apply Dirichlet BCs (`*BOUNDARY`) via elimination/penalty

## 9. Linear solve (spec: linear-algebra-and-solvers)

- [x] 9.1 Solve SPD `K u = f` via SciPP sparse (`spsolve` direct / `cg` SPD-iterative). NOTE: NumPP/SciPP expose no sparse Cholesky; `spsolve`/`cg` are the paths
- [ ] 9.2 Honor `SOLVER=` selection; default solver when unspecified
- [ ] 9.3 Unavailable-solver policy: clear report + documented stop/fallback

## 10. Results (spec: static-analysis, results-output)

- [ ] 10.1 Recover element strains/stresses; extrapolate to nodes; reaction forces
- [ ] 10.2 Write `.frd` (U, S, E, RF) CGX-compatible
- [ ] 10.3 Write `.dat` tabular output honoring output requests

## 11. Python bindings & regression harness (spec: python-bindings, build-and-tooling)

- [ ] 11.1 Bind: build model, load deck, select backend, run, read results as NumPy arrays
- [ ] 11.2 Propagate C++ exceptions as Python exceptions with actionable messages
- [ ] 11.3 Corpus: `beam10p.inp` (C3D10 *CLOAD), `contact4tet.inp` (C3D10 *DLOAD pressure), + hand-authored minimal C3D4 deck; generate reference results with stock CalculiX; per-deck tolerance manifest
- [ ] 11.4 pytest harness: solve each deck, compare U (rel. L2) and S (max rel.) to reference

## 12. CI & docs (spec: build-and-tooling)

- [ ] 12.1 CI: build + `ctest` + pytest regression on the CPU backend
- [ ] 12.2 CI: `openspec validate --all --strict` gate
- [ ] 12.3 Add iOS/Android toolchain files as configure-only CI checks
- [ ] 12.4 Update README/docs with the working end-to-end example
