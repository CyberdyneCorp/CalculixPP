## ADDED Requirements

### Requirement: CMake C++20 build

The project SHALL build with CMake using the C++20 standard and SHALL produce a static/shared solver core library, a command-line runner, and (optionally) the Python module from a single top-level configure step.

#### Scenario: Configure and build the core

- GIVEN a checkout with CMake ≥ 3.24 and a C++20 compiler
- WHEN `cmake -B build && cmake --build build` is run
- THEN the solver core library and the CLI runner SHALL build without error
- AND the C++20 standard SHALL be enforced (no compiler-specific extensions required)

### Requirement: Portable multi-platform targets

The build SHALL support Linux, macOS, and Windows desktop targets, and SHALL support cross-compilation configurations for iOS and Android without source changes to the core.

#### Scenario: Desktop build

- GIVEN a supported desktop toolchain (GCC, Clang, or MSVC)
- WHEN the project is configured and built
- THEN the core library SHALL compile and its unit tests SHALL pass

#### Scenario: Mobile cross-compile configuration

- GIVEN an iOS or Android CMake toolchain file
- WHEN the core library target is configured with that toolchain
- THEN configuration SHALL succeed using only portable C++20 and the declared dependencies

### Requirement: Dependency integration

The build SHALL integrate NumPP (linear algebra/solvers), CyberCadKernel (geometry/meshing), and pybind11 (Python bindings) via pinned versions using submodules or CMake FetchContent, and SHALL fail configuration with a clear message if a required dependency cannot be resolved.

#### Scenario: Dependencies resolved

- GIVEN network/submodule access to the pinned NumPP and CyberCadKernel revisions
- WHEN the project is configured
- THEN NumPP and CyberCadKernel SHALL be made available to the core targets at their pinned versions

#### Scenario: Missing dependency

- GIVEN a required dependency cannot be resolved
- WHEN the project is configured
- THEN configuration SHALL stop and report which dependency and version is missing

### Requirement: Build and run with no GPU toolkit

The default build SHALL require no CUDA, OpenCL, or Metal toolkit, and the full test suite SHALL pass using only the CPU reference backend when no GPU toolkit is installed.

#### Scenario: No GPU toolkit installed

- GIVEN a build environment with no GPU toolkit present
- WHEN the project is built and `ctest` is run
- THEN the build SHALL succeed AND all tests SHALL pass on the CPU reference backend

### Requirement: Unit and regression test harness

The project SHALL provide a GoogleTest/CTest unit harness for C++ components and a pytest harness that drives the solver through the Python bindings against reference decks.

#### Scenario: Run unit tests

- GIVEN a successful build
- WHEN `ctest` is executed
- THEN the registered C++ unit tests SHALL run and report pass/fail per test

#### Scenario: Run regression harness

- GIVEN the Python module is built and reference `test/` decks are available
- WHEN the pytest regression suite is executed
- THEN each deck SHALL be solved and its results compared against the reference within the documented tolerance

### Requirement: Continuous integration gate

CI SHALL, on every pull request, build the core, run the unit and regression suites on the CPU backend, and run `openspec validate --all --strict`, blocking merge on any failure.

#### Scenario: Pull request validation

- GIVEN a pull request against the default branch
- WHEN CI runs
- THEN the build, unit tests, regression suite, and `openspec validate --all --strict` SHALL all pass for the PR to be mergeable
