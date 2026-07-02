# Mesh Processing

## Purpose

Defines geometry import and mesh preparation for CalculixPP through **CyberCadKernel**,
the in-house C++20 OCCT-equivalent CAD kernel. This capability takes CAD geometry
(STEP/IGES/BREP) or an existing analysis mesh and produces a solver-ready finite
element mesh — surface triangulation and tetrahedral volume meshes — which is then
mapped into the mesh-and-model capability (nodes, elements, surfaces). It also owns
geometry healing, mesh quality reporting, and the failure paths for non-manifold or
open geometry.

**Porting Phase:** Phase 1 — Foundation & vertical slice.

## Requirements

### Requirement: CAD geometry import
The system SHALL import CAD geometry in STEP, IGES, and BREP formats through
CyberCadKernel into an in-memory boundary representation (B-rep) of vertices, edges,
faces, and shells. Import SHALL preserve geometric units and report the format,
entity counts, and any unsupported entities encountered.

#### Scenario: Import a STEP solid
- GIVEN a `.step` file describing a single closed solid
- WHEN it is imported through CyberCadKernel
- THEN a B-rep with faces, edges, and vertices SHALL be available for meshing
- AND the imported units and face count SHALL be reported

#### Scenario: Unsupported CAD entity
- GIVEN an IGES file containing an entity type CyberCadKernel cannot represent
- WHEN it is imported
- THEN the entity SHALL be skipped and reported as unsupported without aborting the import

### Requirement: Geometry healing and validation
The system SHALL provide geometry healing that validates a B-rep and repairs common
defects — gaps between faces, sliver faces, and small edges — before meshing. Healing
SHALL run under user-controllable tolerances and SHALL report which defects were found
and repaired. A B-rep that cannot be healed into a valid shell SHALL be reported as a
failure rather than meshed silently.

#### Scenario: Heal a gapped shell
- GIVEN an imported B-rep whose adjacent faces have gaps below the healing tolerance
- WHEN geometry healing runs
- THEN the gaps SHALL be closed and a validation report SHALL list the repairs

#### Scenario: Irreparable geometry
- GIVEN a B-rep with defects that cannot be resolved within the configured tolerance
- WHEN healing runs
- THEN the operation SHALL fail with a report identifying the offending faces/edges

### Requirement: Surface triangulation
The system SHALL tessellate B-rep faces into a triangle surface mesh with
user-controllable size (maximum edge length) and deviation (chordal/angular tolerance)
controls. The resulting triangulation SHALL be watertight across shared edges when the
input B-rep is a valid closed shell.

#### Scenario: Tessellate with deviation control
- GIVEN a healed B-rep and a chordal deviation tolerance
- WHEN surface triangulation runs
- THEN a triangle mesh SHALL be produced whose triangles respect the size and deviation limits
- AND triangles across shared edges SHALL share coincident nodes (watertight)

#### Scenario: Non-closed surface reported
- GIVEN a B-rep or triangulation that is not closed or not manifold
- WHEN triangulation completes
- THEN the mesh SHALL be flagged as open/non-manifold with the free edges identified

### Requirement: Tetrahedral volume meshing
The system SHALL generate a tetrahedral volume mesh — linear `C3D4` or quadratic
`C3D10` elements — from a closed, manifold triangulated surface, under
user-controllable element size and grading controls. Volume meshing SHALL fail with a
diagnostic when the input surface is not a closed manifold.

#### Scenario: Mesh a closed surface to C3D10
- GIVEN a watertight triangulated surface and a target element size, quadratic requested
- WHEN volume meshing runs
- THEN a `C3D10` tetrahedral mesh filling the enclosed volume SHALL be produced

#### Scenario: Open surface cannot be volume meshed
- GIVEN a triangulated surface with free edges (not closed)
- WHEN volume meshing is attempted
- THEN the operation SHALL fail with a diagnostic naming the free edges, and no mesh SHALL be produced

### Requirement: Mesh quality metrics and reporting
The system SHALL compute per-element quality metrics for generated meshes — including
minimum/maximum dihedral angle, aspect ratio, and Jacobian/scaled-Jacobian — and
produce a summary report (min, max, mean, and count of elements below a configurable
quality threshold). Volume meshing SHALL accept a minimum-quality control that rejects
or refines degenerate elements.

#### Scenario: Quality report after meshing
- GIVEN a generated tetrahedral mesh
- WHEN quality metrics are requested
- THEN a report SHALL list per-metric min/max/mean and the number of elements below the threshold

#### Scenario: Degenerate elements flagged
- GIVEN a volume mesh containing tetrahedra below the minimum scaled-Jacobian threshold
- WHEN quality is evaluated
- THEN those elements SHALL be counted and identified in the report

### Requirement: Mapping into the FE model
The system SHALL map a generated mesh into the solver's finite element model consumed
by the **mesh-and-model** capability, emitting unique node numbers with coordinates,
elements with the correct connectivity ordering for their type, and named surfaces
derived from B-rep faces. Node and element numbering SHALL be stable and unique so the
model is directly consumable by loads, boundary conditions, sections, and output
requests.

#### Scenario: Emit nodes, elements, and surfaces
- GIVEN a completed tetrahedral mesh whose source faces carry names
- WHEN it is mapped into the FE model
- THEN unique `*NODE` coordinates and `*ELEMENT` connectivities SHALL be produced for mesh-and-model
- AND named `*SURFACE` groupings SHALL be derived from the originating B-rep faces

### Requirement: Import of existing analysis meshes
The system SHALL import existing CalculiX/Abaqus meshes (nodes and elements from an
input deck) directly into the FE model, bypassing geometry import, healing,
triangulation, and volume meshing entirely. Such meshes SHALL feed the mesh-and-model
capability unchanged.

#### Scenario: Bypass geometry meshing for an existing deck
- GIVEN a CalculiX `*NODE`/`*ELEMENT` deck for an already-meshed model
- WHEN it is loaded
- THEN the nodes and elements SHALL populate mesh-and-model directly without invoking CyberCadKernel meshing

### Requirement: Python-reachable mesh operations
All mesh-processing operations — geometry import, healing, surface triangulation, volume meshing, quality reporting, and mapping into the FE model — SHALL be reachable from the **python-bindings** capability so meshes can be prepared and inspected from Python scripts and tests.

#### Scenario: Prepare a mesh from Python
- GIVEN a STEP file path in a Python script
- WHEN the mesh-processing operations are invoked through the bindings
- THEN import, triangulation, and volume meshing SHALL run and the resulting mesh and quality report SHALL be readable from Python
