## Why

CalculixPP has a complete target-spec baseline but no code. Before any physics
capability can be built, the project needs a portable C++20 skeleton: a build system,
the two heavy dependencies (NumPP, CyberCadKernel) wired in, a CPU compute backend, and
a way to run and validate results from Python. Phase 1 delivers that skeleton **and**
proves it end to end with the simplest complete pipeline — a linear `*STATIC` solve on a
tetrahedral mesh — so every later phase builds on a stack already known to work.

## What Changes

- Introduce the repository build/tooling foundation: CMake (C++20), multi-platform
  targets (Linux, macOS, Windows, iOS, Android), dependency integration for **NumPP**,
  **CyberCadKernel**, and **pybind11**, a GoogleTest/CTest unit harness, a pytest
  regression harness, and a CI gate (build + tests + `openspec validate --all --strict`).
- Implement the **CPU reference compute backend** (via NumPP) and the `ComputeBackend`
  selection surface, with the guarantee that the core builds and passes tests with **no
  GPU toolkit present**.
- Implement the **linear-static vertical slice** across the existing baseline specs:
  parse the Phase-1 keyword subset → build the FE model (tet elements) → assemble
  `K` and `f` → solve `K u = f` with NumPP → recover stresses/strains → write `.frd`/`.dat`.
- Expose the slice through **Python bindings** and validate displacements/stresses
  against the reference CalculiX result for equivalent `test/` decks within a documented
  tolerance.
- No breaking changes (greenfield). GPU backends, nonlinear statics, and other physics
  are explicitly deferred to later phases.

## Capabilities

### New Capabilities

- `build-and-tooling`: the CMake-based build, dependency integration (NumPP,
  CyberCadKernel, pybind11), cross-platform targets, the no-GPU-required build
  guarantee, the unit/regression test harnesses, and the CI validation gate.

### Modified Capabilities

None. Phase 1 **implements** a subset of requirements that already exist in the
target-spec baseline (`input-deck-parsing`, `mesh-and-model`, `element-sections`,
`material-models`, `loads-and-boundary-conditions`, `static-analysis`,
`linear-algebra-and-solvers`, `compute-backend`, `results-output`, `mesh-processing`,
`python-bindings`). Implementing existing requirements is not a spec-level behavior
change, so no delta is created for them; the tasks list references them directly.

## Impact

- **New code:** `src/` C++20 core (model, parser, elements, assembly, solve, results),
  `backend/` (CPU/NumPP), `python/` (pybind11 module), `cmake/`, `tests/` (C++ + pytest),
  `.github/` CI.
- **Dependencies:** NumPP, CyberCadKernel (git submodules or CMake FetchContent);
  pybind11; GoogleTest. No GPU toolkit required.
- **Validation targets:** a curated subset of reference CalculiX `test/*.inp` linear-
  static decks becomes the Phase-1 regression corpus.
- **Docs:** README quickstart (build, run a deck, run Python), and this change archives
  into the baseline once the slice passes.
