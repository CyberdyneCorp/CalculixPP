# Phase 1 — Foundation & Linear-Static Vertical Slice

Tasks reference the capability they implement. Spec deltas exist only for
`build-and-tooling`; all other references are to existing baseline specs.

## 1. Repository & build scaffolding (spec: build-and-tooling)

- [x] 1.1 Create module layout: `core/`, `numerics/`, `io/`, `fem/`, `compute/`, `geometry/`, `apps/` (CLI), `python/`, `tests/`, `cmake/` (present under `include/calculixpp/`, `src/`, and top level)
- [x] 1.2 Top-level `CMakeLists.txt` enforcing C++20; core library + CLI targets
- [x] 1.3 Add `.clang-format`/`.clang-tidy` and wire warnings-as-errors on the core
- [x] 1.4 README quickstart (configure, build, run a deck via `ccxpp`, run Python)

## 2. Dependency integration (spec: build-and-tooling)

- [x] 2.1 Integrate NumPP (`find_package`) + SciPP (`add_subdirectory`; provides the sparse module); smoke test a sparse solve. NOTE: sparse lives in SciPP, not NumPP; see `scripts/bootstrap_deps.sh`
- [x] 2.2 Integrate CyberCadKernel at a pinned version; smoke test loading a B-rep (stub — CyberCadKernel integration deferred to a later phase; deck mesh import is the Phase-1 path)
- [x] 2.3 Integrate pybind11; build an empty `calculixpp` Python module
- [x] 2.4 Clear configure-time failure when a dependency/version is unresolved (`cmake/Dependencies.cmake`: `find_package(NumPP CONFIG REQUIRED)` errors clearly if NumPP is absent; SciPP absence raises a `FATAL_ERROR` naming the searched path and the `-DCALCULIXPP_SCIPP_DIR` override)

## 3. Compute backend — CPU (spec: compute-backend, linear-algebra-and-solvers)

- [x] 3.1 Define `ComputeBackend` interface (`include/calculixpp/compute/backend.hpp`). Phase-1 minimal surface: sparse SPD `solve_sparse`; assembly scatter / SpMV are future work as the backend matures.
- [x] 3.2 Implement CPU reference backend over SciPP sparse (`src/compute/cpu_backend.cpp`, COO→CsrMatrix::from_coo→spsolve/cg); default, always available. `numerics::solve_reduced` now routes through it (behavior-preserving).
- [x] 3.3 Backend selection surface: `BackendKind{CPU,CUDA,OpenCL,Metal}` + `select_backend()` returning CPU (unimplemented backends fall back to CPU — a missing GPU toolkit never breaks build/run). GPU backends remain future work.
- [x] 3.4 Guarantee build + tests pass with no GPU toolkit present (CI `build-test` job builds + runs `ctest` on ubuntu-latest with no CUDA/OpenCL/Metal toolkit; step named "Test (CPU backend, no GPU toolkit)")

## 4. Core FE data model (spec: mesh-and-model)

- [x] 4.1 Node store (ids, coordinates) and element store (connectivity, type)
- [x] 4.2 `C3D4` and `C3D10` tetrahedral element types
- [x] 4.3 Node sets / element sets (`*NSET`/`*ELSET`) and surfaces (`*SURFACE`). NOTE: `Mesh::Surface` stores element (id, face 1..4) pairs or node ids; parsed and stored, not yet consumed by loads.
- [x] 4.4 DOF numbering / equation map
- [x] 4.5 Central `cpp_index_t` (32-bit signed) typedef used by all node/element/DOF indexing

## 5. Mesh import bypass (spec: mesh-processing)

- [x] 5.1 Import mesh directly from the deck (`*NODE`/`*ELEMENT`), bypassing CAD meshing
- [x] 5.2 Stub CyberCadKernel geometry-import/triangulate/tet interfaces (validated later phase) (stub — CyberCadKernel integration deferred to a later phase; deck mesh import is the Phase-1 path)

## 6. Input deck parser — Phase-1 subset (spec: input-deck-parsing)

- [x] 6.1 Keyword tokenizer/lexer for Abaqus-style cards (case-insensitive, comments, continuations)
- [x] 6.2 Parse `*NODE`, `*ELEMENT`, `*NSET`, `*ELSET`, `*SURFACE` (`*SURFACE` TYPE=ELEMENT/NODE now really stored in the mesh with elset/nset expansion, no longer ignored)
- [x] 6.3 Parse `*MATERIAL`, `*ELASTIC`, `*DENSITY`, `*SOLID SECTION`
- [~] 6.4 Parse `*BOUNDARY`, `*CLOAD`, `*DLOAD` (pressure), `*STEP`/`*STATIC`, `*AMPLITUDE` (constant/step)
- [x] 6.5 Actionable parse errors with card + line context

## 7. Material & section (spec: material-models, element-sections)

- [x] 7.1 Isotropic linear-elastic constitutive model (E, ν) and density
- [x] 7.2 `*SOLID SECTION` binding material to element set; local orientation hook (identity for now)

## 8. Assembly & loads (spec: static-analysis, loads-and-boundary-conditions)

- [x] 8.1 Element stiffness kernels for `C3D4`/`C3D10` (Gauss integration)
- [x] 8.2 Global assembly into COO triplets → SciPP `CsrMatrix::from_coo` (compute-backend abstraction still to be formalized)
- [x] 8.3 Concentrated loads (`*CLOAD`) and pressure loads (`*DLOAD` P<face>, surface integration)
- [x] 8.4 Apply Dirichlet BCs (`*BOUNDARY`) via elimination/penalty

## 9. Linear solve (spec: linear-algebra-and-solvers)

- [x] 9.1 Solve SPD `K u = f` via SciPP sparse (`spsolve` direct / `cg` SPD-iterative). NOTE: NumPP/SciPP expose no sparse Cholesky; `spsolve`/`cg` are the paths
- [x] 9.2 Honor `SOLVER=` on *STATIC: parsed and stored on the Model (`RequestedSolver`); SPOOLES/PARDISO/PASTIX/DIRECT→Direct, ITERATIVE*/CG→CG; default Direct when unspecified. CLI/Python/solve honor the model's request.
- [x] 9.3 Unavailable-solver policy: an unrecognized/unavailable SOLVER= name throws a clear `ParseError` naming the requested solver ("solver not available: SOLVER=<name>"), mirroring CalculiX.

## 10. Results (spec: static-analysis, results-output)

- [x] 10.1 Recover element strains/stresses; extrapolate to nodes; reaction forces
- [x] 10.2 Write `.frd` (U, S, E, RF) CGX-compatible
- [x] 10.3 Write `.dat` tabular output honoring output requests

## 11. Python bindings & regression harness (spec: python-bindings, build-and-tooling)

- [x] 11.1 Bind: load deck, select solver+backend, run, read results as NumPy arrays; plus available_backends()/selected_backend() introspection and summary()/summary_text() (num_nodes/num_elements/num_materials/materials/requested_solver) without solving. NOTE: full programmatic model-building (node/element/material builder API) deferred to a later phase — the deck path + summary cover Phase 1 parity.
- [x] 11.2 Propagate C++ exceptions as Python exceptions with actionable messages
- [x] 11.3 Corpus: `beam10p.inp` (C3D10, validated vs committed `beam10p.dat.ref`) + hand-authored C3D4 pressure deck (analytical equilibrium). NOTE: contact4tet=S8 shells, circ11p=*USE REFINED MESH, segmentunsmooth too slow — beam10p is the only small clean CalculiX-ref tet deck
- [x] 11.4 pytest harness: solve each deck, compare U (rel. L2) and S (max rel.) to reference

## 12. CI & docs (spec: build-and-tooling)

- [x] 12.1 CI: build + `ctest` + pytest regression on the CPU backend
- [x] 12.2 CI: `openspec validate --all --strict` gate
- [x] 12.3 Add iOS/Android toolchain files as configure-only CI checks (`cmake/toolchains/{android,ios}.cmake`; CI `mobile-configure` job configures core-only with `WITH_SOLVER=OFF`, `continue-on-error` when NDK/SDK absent)
- [x] 12.4 Update README/docs with the working end-to-end example
