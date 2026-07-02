# Project Context

## Purpose

**CalculiX++ (CalculixPP)** is a ground-up port of the CalculiX CrunchiX (CCX)
three-dimensional structural finite element solver
(https://github.com/Dhondtguido/CalculiX) to **pure, modern C++20**. The goal is a
single, portable, dependency-light solver core that runs on **mobile (iOS, Android)
and desktop**, with **optional hardware acceleration** (CUDA, OpenCL, Metal) behind a
uniform compute abstraction, and **first-class Python bindings** for scripting,
testing, and validation.

The original CalculiX is ~1190 files of mixed Fortran (numerical kernels, element
routines) and C (orchestration), tied to build-time third-party libraries (SPOOLES,
PARDISO, PaStiX, ARPACK, BLAS/LAPACK). CalculixPP replaces that stack with in-house,
portable C++20 components while preserving the physics and Abaqus-style input model
that make CalculiX useful.

### Relationship to the reference implementation

These `openspec/specs/` are **target specifications** for the port. They describe the
behavior CalculixPP SHALL deliver, grounded in the observed behavior of the reference
CalculiX solver. Each requirement is a contract the C++20 implementation must satisfy;
correctness is measured against the reference CalculiX results (the `test/` corpus of
`*.inp` decks) within numerical tolerance. Source breadcrumbs cite the reference
Fortran/C as `(ref: src/<file>)` — these point at the behavior being reproduced, not
code that exists in this repository yet. As capabilities are implemented, changes flow
through the normal OpenSpec change → apply → archive cycle against this baseline.

## Tech Stack

- **Language:** C++20 (concepts, ranges, `std::span`, modules where the toolchains
  allow). No Fortran. No runtime GPU dependency required for a correct build.
- **Linear algebra & solvers:** **NumPP** (https://github.com/CyberdyneCorp/NumPP) —
  in-house C++20 numerics providing dense/sparse containers, factorizations, iterative
  solvers, and eigensolvers. Replaces BLAS/LAPACK, SPOOLES/PARDISO/PaStiX, and ARPACK.
- **Scientific routines:** **SciPP** (https://github.com/CyberdyneCorp/SciPP) — in-house
  C++20 SciPy-equivalent (higher-level numerics on top of NumPP), used where needed for
  quadrature, interpolation, optimization, and sparse/scientific helpers.
- **Geometry & meshing:** **CyberCadKernel**
  (https://github.com/CyberdyneCorp/CyberCadKernel) — a C++20 OCCT-equivalent CAD
  kernel used for geometry import, healing, surface triangulation, and volume
  (tetrahedral) mesh preparation feeding the solver.
- **Compute backends:** a pluggable `ComputeBackend` abstraction. Reference CPU
  backend (via NumPP) is the default and always available; CUDA / OpenCL / Metal
  backends are optional, selected at runtime, and never required for correctness.
- **Python bindings:** **pybind11**, exposing the full public C++ API for scripting,
  automation, and regression testing against reference decks.
- **Build:** CMake (C++20). Platform targets: Linux, macOS, Windows, iOS, Android.
- **Test:** GoogleTest/CTest for C++ units; pytest driving the Python bindings for
  end-to-end validation against the reference `test/` decks.

## Project Conventions

- **Input language:** CalculixPP reads the same Abaqus-style keyword deck as CalculiX
  (`*NODE`, `*ELEMENT`, `*MATERIAL`, `*STATIC`, `*FREQUENCY`, ...). Existing decks are
  the primary compatibility and regression target.
- **Analysis dispatch:** procedures are selected by their keyword card as in CalculiX
  (static, frequency, buckle, modal dynamic, steady-state, heat transfer, CFD/network,
  electromagnetics, contact, crack propagation, sensitivity/optimization).
- **Portability first:** the core solver SHALL compile and run correctly with no GPU
  toolkit present. Acceleration is additive; a backend absence degrades to the CPU
  path, never to an error.
- **Determinism & tolerance:** results SHALL match the reference CalculiX output for
  the corresponding `test/` deck within a documented numerical tolerance.
- **Cognitive complexity:** systems/numerics code target ≤ 25–35 per function
  (element kernels, assembly); bindings and orchestration target ≤ 15. Isolate element
  math behind small, testable kernels.
- **API surface:** every solver capability exposed in C++ SHALL be reachable from the
  Python bindings (full API parity).

## Porting Roadmap (phases)

Specs carry a **Porting Phase** note. Phases gate scope, not spec authorship — all
target specs exist now; implementation lands phase by phase.

- **Phase 1 — Foundation & vertical slice:** build/CMake, NumPP + CyberCadKernel
  integration, compute-backend abstraction (CPU), input-deck parsing subset, tet mesh
  model, linear-elastic solid sections, boundary conditions + concentrated/pressure
  loads, **linear `*STATIC` solve `K u = f`**, `.frd`/`.dat` results, Python bindings.
  This proves the whole pipeline end to end.
- **Phase 2 — Nonlinear statics & materials:** Newton-Raphson, incrementation,
  plasticity, more element families, distributed loads/amplitudes, constraints (MPC,
  equations, ties).
- **Phase 3 — Thermal & coupled:** heat transfer, coupled temperature-displacement,
  contact, model change (element/contact birth-death).
- **Phase 4 — Dynamics & eigenproblems:** modal/buckling, direct & modal dynamics,
  Green-function steps, substructure/superelement generation.
- **Phase 5 — Advanced physics:** CFD/1-D networks, electromagnetics, crack
  propagation, design sensitivity/optimization, adaptive mesh refinement, submodeling,
  high-cycle fatigue.
- **Cross-cutting (per phase):** GPU backends (CUDA/OpenCL/Metal) added where
  assembly/solve hot paths justify it.

## Important Constraints

- **No mandatory third-party numeric libraries.** NumPP and CyberCadKernel are the
  only heavy dependencies, both C++20 and portable.
- **Mobile targets** constrain memory and forbid desktop-only assumptions (no reliance
  on 8-byte-integer builds, no unbounded stack use, no OS-specific solver blobs).
- **Numerical fidelity** to the reference solver is the acceptance criterion for every
  physics capability; refactors must preserve results within tolerance.
- Capabilities that in CalculiX depend on an optionally-linked library map onto NumPP
  functionality; where NumPP lacks a capability, the spec SHALL state the gap rather
  than assume it.

## Domain Context

This is scientific/engineering software whose correctness is defined by physics and
numerical accuracy. Specs describe *what the solver does* (procedures, inputs, outputs,
error conditions), not a re-derivation of the mechanics. The authoritative theory
reference is Dhondt, G., *The Finite Element Method for Three-Dimensional
Thermomechanical Applications*, Wiley, 2004, and the CalculiX user manual.

## External Dependencies

- **NumPP** — C++20 linear algebra, sparse/dense solvers, eigensolvers.
- **SciPP** — C++20 scientific routines (SciPy-equivalent) layered on NumPP.
- **CyberCadKernel** — C++20 CAD/geometry kernel for import and mesh preparation.
- **pybind11** — Python binding layer.
- **CUDA / OpenCL / Metal** — optional acceleration backends (runtime-selected).
- Reference only (not linked): CalculiX (Fortran/C) as the behavioral oracle and its
  `test/` deck corpus; CGX `.frd` format as the results interchange target.
