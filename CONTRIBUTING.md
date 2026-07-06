# Contributing to CalculiX++

Thanks for your interest in CalculiX++ (`CalculixPP`) — a ground-up port of the CalculiX 3-D
structural finite element solver to modern, portable C++20. This guide covers how to build, the
workflow we follow, and what a mergeable change looks like.

## Getting set up

```bash
git clone https://github.com/CyberdyneCorp/CalculixPP.git
cd CalculixPP
just bootstrap      # build + install NumPP into .deps/ and point at a SciPP checkout
just test           # configure + build + full CTest suite (CPU backend)
```

No `just`? The equivalent is `scripts/bootstrap_deps.sh --numpp ../NumPP --scipp ../SciPP` then:

```bash
cmake -S . -B build -G Ninja \
  -DCALCULIXPP_WITH_SOLVER=ON -DCALCULIXPP_BUILD_PYTHON=ON \
  -DCMAKE_PREFIX_PATH=$PWD/.deps/install -DCALCULIXPP_SCIPP_DIR=../SciPP
cmake --build build
ctest --test-dir build --output-on-failure
```

Requires a C++20 compiler and CMake ≥ 3.24. The dependency-free `calculixpp_core` also builds
standalone with `-DCALCULIXPP_WITH_SOLVER=OFF` (no NumPP/SciPP) — that is the path the mobile
cross-compile toolchains use. **No GPU toolkit is required** (`just gpu-detect` reports what your
host has). Python bindings additionally need `pip install pybind11 pytest numpy`.

## Development workflow

We develop spec-first with [OpenSpec](https://openspec.dev). Living capability specs are in
[`openspec/specs/`](openspec/specs); completed changes are archived under
[`openspec/changes/archive/`](openspec/changes/archive).

- **Medium or large features** (new physics, elements, materials, procedures, backends, bindings,
  build/packaging behavior): start with an OpenSpec change (proposal → specs/tasks) before
  implementing. Run `openspec validate --all --strict` — CI enforces it.
- **Small fixes and docs:** a direct PR is fine. Still keep the specs and docs truthful — if your
  change makes a spec or doc statement false, update it in the same PR.

## What a mergeable PR looks like

- **Tests.** New behavior ships with tests. **Every bug fix includes a regression test** that
  fails before the fix and passes after.
- **Green CI.** The OpenSpec validation gate, the Linux build + CTest suite (core, numerics, and
  Python regression), and the mobile cross-compile configure smoke check must all pass.
- **Docs/specs in sync.** Update the `README.md` and `openspec/specs/` when your change affects
  documented behavior. Don't leave a claim the code contradicts.
- **Readable, low-complexity code.** Match the surrounding style (C++20 — Concepts, Ranges,
  `std::span`, smart pointers, contiguous data; no raw owning pointers, no legacy macros). Keep
  per-function cognitive complexity modest; isolate genuinely irreducible element/numerics kernels
  and flag them rather than mangling them to hit a number.
- **A descriptive PR message.** Explain what changed and why. If you found a bug, describe how it
  reproduced.

## Commit & PR style

- Concise, technical, imperative commit subjects (e.g. "Fix C3D20 shape-function Jacobian sign").
- Reference the OpenSpec change or issue when there is one.
- Keep unrelated changes in separate PRs.

## Numerical correctness

Physics changes are validated against **stock CalculiX** as the behavioral oracle (the `test/`
deck corpus — e.g. `beam10p`, `beam8p`, `beam20p`, `beamt`, `beam8f`, `beamb`) or against an
analytic reference, within a documented numerical tolerance. If you change an element kernel, a
material model, a constraint, or a solve path, make sure the relevant reference-deck / analytic
tests still pass — and add one if a gap exists.

## Reporting bugs & requesting features

Open an issue using the templates. For **security** issues, do **not** open a public issue — see
[SECURITY.md](SECURITY.md).

## License

By contributing, you agree that your contributions are licensed under the project's
[MIT License](LICENSE).
