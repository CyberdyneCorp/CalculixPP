# Changelog

All notable changes to CalculiX++ are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project follows
[Semantic Versioning](https://semver.org/spec/v2.0.0.html) (pre-1.0: minor/patch may still carry
breaking changes while the API stabilizes).

## [Unreleased]

## [0.2.0] - 2026-07-06

### Added
- **Installable package** — `install`/`export` rules and a generated `CalculixPPConfig.cmake` so
  downstream projects can `find_package(CalculixPP)` and link `CalculiXPP::core`, the
  dependency-free solver core (domain model + FE element math + assembly + `.inp`/`.frd` I/O).
  Ships a `CalculixPPConfigVersion.cmake` with `SameMajorVersion` compatibility.
- **`justfile`** — developer task runner (`build`, `test`, `python`, `core`, `debug`, `spec`,
  `ci`, `clean`) and a cross-platform `gpu-detect` recipe (CUDA / OpenCL / Metal probe).
- **Project governance** — `LICENSE` (MIT), `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`,
  `SECURITY.md` (private vulnerability reporting), this changelog, and GitHub issue/PR templates.
- **Consumability spec** — a portable OpenSpec "usable by others" readiness rubric
  (`openspec/specs/consumability`) consolidating the install / packaging / CI / versioning /
  governance requirements, reusable as an adoption checklist for other projects.
- **README** — a "Versioning & API stability" section, a Contributing section, and installed-core
  `find_package(CalculixPP)` usage.

### Changed
- `calculixpp_core` include directories now use `BUILD_INTERFACE`/`INSTALL_INTERFACE` generator
  expressions so the target can be exported for `find_package`.

## [0.1.0] - 2026-07-05

Initial development version — a ground-up C++20 port of the CalculiX 3-D structural finite element
solver, validated against stock CalculiX references.

### Added
- **Phase 1 — Foundation** — build system, NumPP/SciPP integration, CPU compute backend,
  Abaqus-style `.inp` parsing, tet/hex element math, linear-elastic `*STATIC` solve (`K u = f`,
  sparse direct / IC0-CG), CGX-compatible `.frd`/`.dat` output, and pybind11 Python bindings.
- **Phase 2 — Nonlinear statics & materials** — Newton-Raphson driver with automatic
  incrementation, J2 plasticity, neo-Hookean hyperelasticity, C++ user material, the hex/wedge
  element family, amplitudes/body loads, and constraints (`*EQUATION`/`*MPC`/`*RIGID BODY`/
  `*COUPLING`/`*TIE`).
- **Phase 3 — Thermal & contact** — steady/transient heat transfer, coupled
  temperature-displacement, cavity radiation, node-to-surface penalty contact with a spatial
  contact-search engine, Coulomb friction, thermal contact, and multi-step model change.
- **Phase 4 — Dynamics & eigenproblems** — mechanical mass matrix, generalized eigensolution
  (`*FREQUENCY`), modal/direct/steady-state dynamics, Craig-Bampton/Guyan substructures,
  `*BUCKLE`, and proportional-damping `*COMPLEX FREQUENCY`.
- **Phase 5 — Structural completion (unblocked)** — design sensitivity/optimization (adjoint),
  submodeling, and high-cycle fatigue.
- **Tooling** — a GitHub Actions CI pipeline (OpenSpec validation + Linux build/tests + mobile
  cross-compile configure) and in-repo OpenSpec specs.

[Unreleased]: https://github.com/CyberdyneCorp/CalculixPP/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/CyberdyneCorp/CalculixPP/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/CyberdyneCorp/CalculixPP/releases/tag/v0.1.0
