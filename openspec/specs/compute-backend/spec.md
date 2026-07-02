# Compute Backend

## Purpose

Defines a pluggable `ComputeBackend` abstraction that hosts the numerical hot paths
(assembly, sparse matrix-vector products, factorization) behind a uniform interface. A
CPU reference backend built on NumPP is always available and is the default; optional
CUDA, OpenCL, and Metal backends are runtime-selected and never required for a correct
build or run. Mobile-first portability is paramount: absence of any GPU toolkit degrades
to the CPU path, never to an error. Consumed by linear-algebra-and-solvers and the
assembly capabilities.

**Porting Phase:** Phase 1 (CPU reference backend); GPU backends are cross-cutting and added per later phases where hot paths justify it.

## Requirements

### Requirement: Pluggable backend with CPU default
CalculixPP SHALL expose a `ComputeBackend` abstraction whose CPU reference backend, built
on NumPP, is always available and is the default. Optional CUDA, OpenCL, and Metal
backends SHALL implement the same interface and be selectable at runtime, discovered and
queried by capability, so a consumer requests an operation without knowing the concrete
backend.

#### Scenario: CPU backend is the default
- GIVEN no backend has been explicitly selected
- WHEN a compute capability (assembly, spmv, or factorization) is requested
- THEN the CPU reference backend via NumPP SHALL be used

#### Scenario: Capability query
- GIVEN a runtime-selected optional backend
- WHEN a consumer queries whether it supports a given capability
- THEN the backend SHALL report its supported capabilities so the consumer can offload only what is available

### Requirement: Build and test with no GPU toolkit
The core SHALL compile and pass the full test suite with no CUDA, OpenCL, or Metal
toolkit present. GPU support is additive; its absence SHALL NOT break the build or a
correct run on iOS, Android, or desktop.

#### Scenario: No GPU toolkit installed
- GIVEN a build environment with no GPU toolkit installed
- WHEN the project is built and the test suite is run
- THEN the build SHALL succeed AND all tests SHALL pass using the CPU reference backend

### Requirement: CPU multithreading with configurable thread count
The CPU reference backend SHALL execute matrix assembly and the solve multithreaded, and the thread count SHALL be configurable at runtime — via an explicit API setting and via an `OMP_NUM_THREADS`-style environment variable — defaulting to the hardware concurrency when unset. This replaces CalculiX's OpenMP thread control (`OMP_NUM_THREADS`, `CCX_NPROC_*`). (ref: src/CalculiX.c, src/pthread wrappers)

#### Scenario: Explicit thread count honored
- GIVEN a requested thread count set via the API or the environment variable
- WHEN assembly and the solve run
- THEN the CPU backend SHALL use the requested number of threads

#### Scenario: Default to hardware concurrency
- GIVEN no thread count is configured
- WHEN a compute operation runs
- THEN the CPU backend SHALL default to the available hardware concurrency

#### Scenario: Thread count does not change results
- GIVEN the same system solved with different thread counts
- WHEN the solves complete
- THEN the results SHALL agree within the documented numerical tolerance

### Requirement: Backend selection policy with graceful fallback
Backend selection SHALL follow a documented order: an explicit user choice first, else
the best available accelerated backend, else the CPU reference backend. When an explicitly
requested or auto-selected device is missing or unsupported, selection SHALL fall back
gracefully to the next option and ultimately to CPU, emitting a diagnostic rather than
failing the run.

#### Scenario: Explicit choice honored
- GIVEN a user explicitly selects an available backend
- WHEN a compute capability runs
- THEN that backend SHALL be used

#### Scenario: Graceful fallback when device is missing
- GIVEN a user requests a backend whose device is absent or unsupported at runtime
- WHEN selection resolves
- THEN CalculixPP SHALL emit a diagnostic AND fall back to the next available backend, ending at the CPU reference backend, so the run still completes correctly

### Requirement: Results match the CPU reference within tolerance
Numerical results produced by any backend SHALL match the CPU reference backend within a
documented tolerance. A backend that cannot meet the tolerance for a capability SHALL not
be selected for that capability.

#### Scenario: Accelerated result matches CPU
- GIVEN the same input solved by an accelerated backend and by the CPU reference backend
- WHEN both complete
- THEN the results SHALL agree within the documented tolerance

### Requirement: Defined memory and host-device transfer boundaries
Device memory management and host-device transfer boundaries for assembled matrices and vectors SHALL be defined so that a capability (assembly, spmv, or factorization) can offload to a device without leaking backend-specific details into the caller. Ownership, lifetime, and transfer points for the shared sparse data SHALL be explicit.

#### Scenario: Offload without leaking backend details
- GIVEN an assembled sparse matrix and vectors in the shared host representation
- WHEN a capability offloads to a device backend
- THEN the backend SHALL manage device allocation and host-device transfers internally, and the caller SHALL interact only through the abstraction

#### Scenario: Deterministic release of device memory
- GIVEN a capability that allocated device memory for a solve
- WHEN the capability completes or is torn down
- THEN the associated device memory SHALL be released without the caller managing device handles

### Requirement: Backend selection from Python bindings
Backend discovery and selection SHALL be reachable from the Python bindings, at full
parity with the C++ API, so scripts and regression tests can enumerate available backends
and pin a specific one.

#### Scenario: Select backend from Python
- GIVEN a Python script using the pybind11 bindings
- WHEN it enumerates available backends and selects one (or requests the default)
- THEN the chosen backend SHALL be used for subsequent solves, with the same fallback policy as the C++ API
