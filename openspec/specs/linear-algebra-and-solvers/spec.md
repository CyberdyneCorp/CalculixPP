# Linear Algebra and Solvers

## Purpose

Solves the sparse linear systems assembled by every procedure (`K u = f` and its
variants) and the generalized eigenproblems used by modal analyses. In CalculixPP the
solver layer is a thin CalculiX-facing abstraction over **NumPP**, replacing the
reference stack of SPOOLES / PARDISO / PaStiX (direct), the CalculiX iterative solver
(preconditioned CG), and ARPACK (eigenproblems) with in-house, portable C++20 numerics.
All solves execute through the **ComputeBackend** so acceleration is additive.
(ref: src/spooles.c, src/pardiso.c, src/pastix.c, src/preiter.c, src/arpack.c)

**Porting Phase:** Phase 1 (direct/iterative sparse solves for linear `*STATIC`); eigensolvers land with Phase 4 (see compute-backend, modal-and-buckling-analysis).

## Requirements

### Requirement: Selectable solver over NumPP
The solver interface SHALL be a CalculixPP abstraction over NumPP, selected per step via
the `SOLVER` parameter on the procedure card (e.g. `*STATIC, SOLVER=SPOOLES|PARDISO|
PASTIX|ITERATIVE SCALING|ITERATIVE CHOLESKY`). When `SOLVER` is unspecified, a NumPP
direct sparse solver SHALL be used as the default. Legacy solver names SHALL be accepted
as aliases mapping to the equivalent NumPP capability. (ref: src/calinput.f, src/CalculiX.c)

#### Scenario: Default solver when unspecified
- GIVEN a `*STATIC` step with no `SOLVER` parameter
- WHEN the assembled system `K u = f` is solved
- THEN CalculixPP SHALL use the NumPP direct sparse solver as the default

#### Scenario: Explicit legacy solver name
- GIVEN a step requesting `SOLVER=SPOOLES` or `SOLVER=PARDISO`
- WHEN the step runs
- THEN the name SHALL be mapped to the equivalent NumPP direct factorization capability

### Requirement: Direct sparse factorization
CalculixPP SHALL provide direct sparse factorization via NumPP for both symmetric and
non-symmetric systems, obtaining the solution by back-substitution. This replaces the
reference SPOOLES / PARDISO / PaStiX direct solvers. (ref: src/spooles.c, src/pardiso.c, src/pastix.c)

#### Scenario: Symmetric direct factorization
- GIVEN a step producing a symmetric system and a direct solver selection
- WHEN the linear system is solved
- THEN NumPP SHALL factorize the matrix directly and solve by back-substitution

#### Scenario: Non-symmetric direct factorization
- GIVEN a procedure that produces a non-symmetric system
- WHEN a direct solver is selected
- THEN NumPP SHALL use a factorization appropriate to non-symmetric matrices and solve to convergence

### Requirement: Iterative solvers with preconditioning
CalculixPP SHALL provide preconditioned iterative solvers via NumPP — conjugate gradient
for symmetric positive-definite systems and GMRES for general systems — with scaling and
incomplete-Cholesky / ILU preconditioners, for large sparse systems. (ref: src/preiter.c, src/iter.f)

#### Scenario: Preconditioned CG on an SPD system
- GIVEN a large symmetric positive-definite system and `SOLVER=ITERATIVE CHOLESKY` (or `ITERATIVE SCALING`)
- WHEN the system is solved
- THEN NumPP SHALL run a preconditioned conjugate-gradient iteration (scaling or incomplete-Cholesky preconditioner) instead of a direct factorization

#### Scenario: GMRES on a general system
- GIVEN a large non-symmetric system and an iterative solver selection
- WHEN the system is solved
- THEN NumPP SHALL run a preconditioned GMRES iteration with an ILU preconditioner

### Requirement: Generalized symmetric eigensolver
CalculixPP SHALL solve the generalized symmetric eigenproblem `K x = λ M x` via NumPP,
replacing ARPACK, for the frequency, buckling, and complex-frequency procedures. Eigen
extraction is cross-referenced by the modal-and-buckling-analysis capability. (ref: src/arpack.c, src/arpackbu.c, src/arpackcs.c)

#### Scenario: Frequency / buckling eigenpairs
- GIVEN a `*FREQUENCY` or `*BUCKLE` step requesting a number of modes
- WHEN eigenpairs are requested
- THEN NumPP SHALL compute the requested eigenvalues and eigenvectors of the generalized symmetric problem

### Requirement: Dense linear algebra for element-level problems
Small dense problems — element stiffness/mass formation, principal-stress and small symmetric eigenproblems (e.g. `DSYEV`-class), and small dense linear solves (`DGESV`/`DGELSS`-class) — SHALL be performed through NumPP dense routines, which SHALL be available in every build (replacing the LAPACK/BLAS routines linked unconditionally by reference CalculiX). (ref: src/calceigenvalues.f, LAPACK DSYEV/DGESV/DGELSS, reference BLAS)

#### Scenario: Small symmetric eigenproblem for principal stresses
- GIVEN a symmetric stress tensor at an integration point
- WHEN principal stresses/directions are computed
- THEN CalculixPP SHALL obtain them via a NumPP dense symmetric eigensolve available in every build

#### Scenario: Element-level dense solve
- GIVEN an element routine requiring a small dense linear solve or least-squares fit
- WHEN the element quantity is formed
- THEN the solve SHALL use NumPP dense routines rather than a hand-rolled factorization

### Requirement: Shared sparse storage format
Assembled matrices SHALL use a compressed sparse (CSR-style) storage format shared
between the assembly step and the NumPP solver interface, so that assembly, factorization,
and sparse matrix-vector products operate on a single representation without conversion
copies beyond what a backend transfer requires.

#### Scenario: Assembly to solve without reformatting
- GIVEN a stiffness matrix assembled into the shared compressed-sparse format
- WHEN it is handed to a direct or iterative solver
- THEN the solver SHALL consume the same CSR-style representation without an intermediate reformat step

### Requirement: Unavailable solver policy
When a requested solver or preconditioner is not available in the current NumPP build, CalculixPP SHALL report the unavailable selection clearly and either stop, or fall back according to a documented policy — mirroring the reference behavior where an unlinked solver produces an error. (ref: src/CalculiX.c)

#### Scenario: Requested solver unavailable
- GIVEN a step requests a solver or preconditioner not provided by the current NumPP build
- WHEN the step runs
- THEN CalculixPP SHALL emit a clear diagnostic naming the unavailable selection AND SHALL either stop or apply the documented fallback (never silently produce wrong results)

### Requirement: Solves execute through the ComputeBackend
All factorizations, back-substitutions, and sparse matrix-vector products SHALL be
dispatched through the ComputeBackend abstraction (see compute-backend), so that a GPU
backend can accelerate a solve without changing results beyond the documented tolerance
and without a GPU toolkit being required for a correct run.

#### Scenario: CPU and accelerated solve agree
- GIVEN the same assembled system solved on the CPU reference backend and on an available GPU backend
- WHEN both solves complete
- THEN the two solutions SHALL agree within the documented numerical tolerance
