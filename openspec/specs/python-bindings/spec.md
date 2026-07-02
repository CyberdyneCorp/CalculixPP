# Python Bindings

## Purpose

Defines the Python binding layer for CalculixPP, built with **pybind11**, that exposes
the full public C++ solver API to Python for scripting, automation, and regression
testing. The bindings are the project's primary validation harness: they drive the
reference `test/` decks end to end and compare results against reference CalculiX
within tolerance. Python is a desktop scripting/test layer with full API parity; the
C++ core remains mobile-capable and the Python module is not built for mobile targets.

**Porting Phase:** Phase 1 — Foundation & vertical slice.

## Requirements

### Requirement: Full public C++ API exposure
The Python module SHALL expose the full public C++ surface so a user can, entirely
from Python, build a model programmatically (nodes, elements, sets, surfaces), load an
Abaqus-style keyword deck, define materials/sections/loads/boundary conditions, select
the analysis procedure, choose the compute backend, run the solve, and read results.

#### Scenario: End-to-end model from Python
- GIVEN a Python session importing the CalculixPP module
- WHEN a model is built or a deck is loaded, a procedure is selected, and the solve is run
- THEN the solve SHALL complete and results SHALL be retrievable through the Python API

#### Scenario: Programmatic model construction
- GIVEN Python calls that add nodes, elements, a material, a section, a load, and a boundary condition
- WHEN the model is assembled and solved
- THEN the result SHALL match an equivalent model loaded from a keyword deck within tolerance

### Requirement: NumPy result and mesh arrays
Result fields and mesh arrays (node coordinates, element connectivity, displacements, stresses, reaction forces) SHALL be exposed as NumPy arrays, zero-copy where it is safe to alias the C++ buffers and copying otherwise. The exposed arrays SHALL carry correct shapes and dtypes for direct use in validation and scripting.

#### Scenario: Read displacements as a NumPy array
- GIVEN a completed static solve
- WHEN nodal displacements are requested from Python
- THEN a NumPy array of shape `(n_nodes, n_dof)` SHALL be returned for numerical comparison

#### Scenario: Zero-copy where safe
- GIVEN a large read-only result buffer owned by the C++ core
- WHEN it is exposed to Python
- THEN the NumPy array SHALL alias the C++ memory without copying, and its lifetime SHALL remain valid while referenced

### Requirement: Exception propagation
C++ exceptions and error conditions SHALL propagate to Python as Python exceptions
carrying actionable messages that identify the offending input or operation. A failure
in the C++ core SHALL NOT crash the interpreter.

#### Scenario: Invalid input raises in Python
- GIVEN a deck or API call with an invalid element connectivity
- WHEN it is processed through the bindings
- THEN a Python exception SHALL be raised with a message identifying the offending element, and the interpreter SHALL stay alive

### Requirement: Reference deck regression harness
The bindings SHALL drive the reference `test/` decks end to end from Python (pytest)
and compare computed nodal/element results against the reference CalculiX results for
each deck within a documented numerical tolerance. This harness is the primary
regression mechanism for the port.

#### Scenario: Deck matches reference within tolerance
- GIVEN a reference `test/` deck and its known-good CalculiX results
- WHEN the deck is solved through the Python bindings
- THEN the computed results SHALL match the reference within the documented tolerance, else the test SHALL fail

### Requirement: Compute-backend selection from Python
Compute-backend selection and query SHALL be available from Python, allowing a script
to list available backends, choose one for a solve, and read back the active backend.
Backend selection SHALL follow the **compute-backend** capability's rules — the CPU
backend is always available and an unavailable backend SHALL degrade to CPU, never
error.

#### Scenario: Select and query a backend
- GIVEN a Python script that lists backends and selects one before solving
- WHEN the solve runs
- THEN the active backend SHALL be queryable from Python and SHALL reflect the selection (or CPU fallback when the requested backend is unavailable)

### Requirement: API parity enforcement
API parity SHALL be enforced: every public C++ capability SHALL have a corresponding
Python entry point, and any public C++ capability without one SHALL be treated as a
defect. Parity SHALL be checked by an automated API-coverage test that enumerates the
public C++ surface and asserts each element is reachable from Python.

#### Scenario: Missing binding is a test failure
- GIVEN a public C++ method with no exposed Python entry point
- WHEN the API-coverage test runs
- THEN the test SHALL fail and identify the unbound C++ symbol

### Requirement: Desktop build scope
The Python module SHALL build for desktop platforms (Linux, macOS, Windows) via CMake
and pybind11. Mobile embedding of the Python layer is out of scope; the C++ core SHALL
remain buildable and correct for mobile targets independently of the Python bindings.

#### Scenario: Desktop build produces an importable module
- GIVEN a desktop CMake build with Python bindings enabled
- WHEN the build completes
- THEN an importable Python module SHALL be produced, while a mobile core build SHALL succeed with the Python layer excluded
