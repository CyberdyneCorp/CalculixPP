## Context

CalculixPP is a greenfield C++20 port of CalculiX with a complete target-spec baseline
but no implementation. The reference solver is Fortran/C tied to build-time libraries
(SPOOLES/PARDISO/PaStiX/ARPACK/BLAS/LAPACK); CalculixPP replaces those with NumPP
(numerics) and CyberCadKernel (geometry/meshing), targets mobile as well as desktop, and
must remain correct with no GPU toolkit present. This change establishes the skeleton and
proves it with the linear-static slice. Correctness is defined against reference CalculiX
`test/` decks within a documented numerical tolerance.

## Goals / Non-Goals

**Goals:**

- A CMake C++20 build that compiles the core, CLI runner, and Python module.
- NumPP + CyberCadKernel + pybind11 integrated at pinned versions.
- A `ComputeBackend` seam with a CPU reference backend; core builds/passes tests with no GPU.
- End-to-end linear `*STATIC` on tets: parse subset → assemble → NumPP solve → stresses → `.frd`/`.dat`.
- Python bindings for the slice and a pytest regression harness comparing to reference results.
- CI gate: build + tests + `openspec validate --all --strict`.

**Non-Goals:**

- GPU backends (CUDA/OpenCL/Metal) — seam only, no implementation this phase.
- Nonlinear statics, plasticity, contact, dynamics, thermal, CFD, EM, optimization.
- Full element library, full keyword coverage, full material library.
- Full Python API parity (only the slice's surface is bound now; parity is enforced later).

## Decisions

- **Module layout.** `core` (model, parser, elements, assembly, results), `numerics`
  (ComputeBackend + NumPP adapters), `geometry` (CyberCadKernel adapters), `cli`,
  `python`. Rationale: keeps element math and numerics behind narrow interfaces (aligns
  with the cognitive-complexity targets and the ComputeBackend abstraction). Alternative
  considered: monolithic `src/` mirroring CalculiX — rejected as it re-imports the
  original's coupling.

- **ComputeBackend seam now, GPU later.** Define `ComputeBackend` with the operations the
  slice needs (assembly scatter, SpMV, sparse factor/solve) and implement only the CPU
  path over NumPP. Rationale: locking the seam early prevents a GPU retrofit from touching
  every call site. Alternative: call NumPP directly and add the seam later — rejected,
  it would require a wide refactor across assembly/solve.

- **NumPP for all algebra.** Sparse storage (CSR), the direct sparse solver (default), and
  the CG/GMRES iterative path all come from NumPP; no BLAS/LAPACK/SPOOLES. Rationale: the
  portability and single-dependency mandate. Trade-off: NumPP solver maturity gates the
  slice — mitigated by starting with the direct solver on small SPD systems.

- **CyberCadKernel for geometry; deck-mesh bypass for the slice.** The slice reads a mesh
  directly from the Abaqus deck (`*NODE`/`*ELEMENT`), exercising the
  `mesh-processing` "existing-mesh import bypass" path; CAD import→triangulate→tet meshing
  is wired as an interface but its full pipeline is validated in a later phase. Rationale:
  fastest route to an end-to-end solve without blocking on the mesher.

- **Elements.** Implement `C3D4` (linear tet) and `C3D10` (quadratic tet) with isotropic
  linear elasticity and a `*SOLID SECTION`. Rationale: tets are the target of the
  mesh-processing pipeline and the simplest complete structural kernel.

- **Results format.** Emit CGX-compatible `.frd` (U, S, E, RF) and `.dat`. Rationale: reuse
  the reference postprocessing ecosystem and enable direct numerical diffing.

- **Tolerance definition.** Agreement is measured as relative L2 norm on nodal
  displacements and max relative error on von Mises stress at integration points, with the
  threshold recorded per deck in the regression manifest. Rationale: displacement L2 is
  stable; stress needs a looser, explicit bound. Alternative: exact float match — rejected
  (different solver/order-of-operations).

- **Python via pybind11.** Bind model build, deck load, backend select, run, and result
  access as NumPy arrays. Rationale: pytest is the regression harness and the fastest
  validation loop.

## Risks / Trade-offs

- **NumPP capability gaps** (missing solver/preconditioner) → start with the direct sparse
  solver; the `linear-algebra-and-solvers` unavailable-solver policy defines the failure
  path; log gaps as follow-up changes.
- **CyberCadKernel/NumPP build friction across mobile toolchains** → validate desktop
  first; add iOS/Android toolchain files as configure-only checks in CI, defer on-device
  runs.
- **Reference-tolerance flakiness** → per-deck thresholds in a manifest, not a single
  global constant; quarantine decks that need investigation rather than loosening globally.
- **Scope creep into Phase 2** → the spec delta is limited to `build-and-tooling`; physics
  tasks reference existing baseline requirements only, so "just add nonlinear" would
  require a new change.

## Migration Plan

Greenfield — no migration. Rollout: land build/tooling, then the slice behind CI; the
change archives into the baseline once the regression corpus passes. Rollback is deleting
the branch (no consumers yet).

## Resolved Decisions (formerly open questions)

- **Dependency delivery → CMake FetchContent with pinned git tags.** NumPP and
  CyberCadKernel are fetched by pinned tag so a single `cmake -B build` configures the
  whole tree with no `submodule init` step — important for the mobile cross-compile
  toolchains and CI. A `-DCPP_USE_SYSTEM_DEPS=ON` (find_package) escape hatch supports
  vendored/offline builds. Rejected plain submodules (extra manual step, easy to forget)
  and unpinned FetchContent (non-reproducible).

- **Minimal Phase-1 reference corpus → three decks.**
  - `beam10p.inp` — C3D10 cantilever, `*CLOAD`; the primary clean end-to-end deck (uses
    only the Phase-1 keyword subset, no out-of-scope cards).
  - `contact4tet.inp` — C3D10, `*DLOAD` pressure; validates distributed loads. (No
    `*CONTACT` card despite the file name — verified; it is a pure static pressure case.)
  - A **hand-authored minimal C3D4 deck** (few linear tets, `*CLOAD` + a Dirichlet BC):
    the reference `test/` corpus has no small pure-linear-tet static deck, so we author
    one and generate its reference result with stock CalculiX.
  Decks using out-of-scope cards (`*REFINE MESH`, contact, `NLGEOM`, plasticity) are
  excluded this phase. Each deck's tolerance lives in the regression manifest (task 11.3).

- **Solver scope → symmetric positive-definite only.** Linear elasticity with adequate
  Dirichlet BCs yields an SPD `K`, so Phase 1 uses NumPP's symmetric sparse direct
  factorization (Cholesky/LDLᵀ). Non-symmetric and indefinite systems are out of scope;
  a request that would produce one is reported via the `linear-algebra-and-solvers`
  unavailable/unsupported-solver path rather than silently mis-solved.

- **Index width → 32-bit signed, via a central `cpp_index_t` typedef.** One 32-bit index
  type for nodes/elements/DOF (no CalculiX-style i8 build); it covers mobile-scale models
  (< 2.1B DOF) and halves index memory versus 64-bit. Because every consumer uses the
  central typedef, widening to 64-bit later is a one-line change if a target model ever
  needs it.
