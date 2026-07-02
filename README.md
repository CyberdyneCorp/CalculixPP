# CalculiX++ (CalculixPP)

A ground-up port of the [CalculiX](https://github.com/Dhondtguido/CalculiX)
three-dimensional structural finite element solver to **modern, portable C++20**,
targeting mobile (iOS/Android) and desktop with optional GPU acceleration.

Numerics come from in-house libraries — **[NumPP](https://github.com/CyberdyneCorp/NumPP)**
(NumPy-equivalent arrays + dense linear algebra) and
**[SciPP](https://github.com/CyberdyneCorp/SciPP)** (SciPy-equivalent; sparse matrices
and sparse solvers) — with **[CyberCadKernel](https://github.com/CyberdyneCorp/CyberCadKernel)**
for CAD/meshing (later phases). The specification and phased roadmap live in
[`openspec/`](openspec/).

## Status

**Phase 1 (Foundation) — in progress.** The dependency-free solver core is up:
domain model, tetrahedral element kernels (C3D4/C3D10 linear elasticity), and global
assembly with single-point-constraint elimination — all unit-tested. Next: the SciPP
sparse solve, the `.inp` parser, stress recovery, `.frd`/`.dat` output, and Python
bindings. See `openspec/changes/phase-1-foundation/tasks.md`.

## Architecture

- **`calculixpp_core`** — dependency-free: domain model + FE element math + assembly to
  COO triplets. Builds and tests everywhere (including mobile toolchains).
- **numerics layer** (next) — hands the assembled triplets to SciPP
  (`scipp::sparse::CsrMatrix` → `spsolve`/`cg`) for the `K u = f` solve.

## Build

```bash
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build --output-on-failure
```

The core requires only a C++20 compiler — **no GPU toolkit, no external libraries**.

### Dependencies (numerics layer, upcoming)

NumPP is consumed via `find_package(NumPP)`; SciPP via `add_subdirectory` (it has no
install/export config) and pulls NumPP transitively. `scripts/bootstrap_deps.sh` builds
and installs NumPP into a local prefix and points the build at a SciPP checkout.

## License

See [LICENSE](LICENSE).
