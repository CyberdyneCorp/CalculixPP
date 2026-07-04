# Mesh and Model Definition

## Purpose

Defines the in-memory finite element model CalculixPP builds from a parsed deck: nodes,
elements (the element library), and the named groupings (sets and surfaces) used to
apply data and requests to portions of the model. Geometry import and volume-mesh
generation (surface triangulation, tetrahedral meshing) are handled upstream by
**CyberCadKernel** and covered by the mesh-processing spec; this capability covers the
resulting FE model once meshed. (ref: src/calinput.f, src/nodes.f, src/elements.f,
src/noset.f, src/elset.f, src/surfaces.f)

**Porting Phase:** 1 — Foundation & vertical slice (nodes, tetrahedral elements C3D4/C3D10, `*NSET`/`*ELSET`, `*SURFACE`; the wider element library and 1-D/2-D expansion land in later phases).

## Requirements

### Requirement: Node definition
The model SHALL define nodes via `*NODE`, each with an integer node number and up to
three Cartesian coordinates. A node number SHALL be unique; redefining it overwrites
the coordinates. Node storage SHALL use portable fixed-width integer indexing suitable
for mobile targets (no 8-byte-integer build assumption). (ref: src/nodes.f)

#### Scenario: Define a node
- GIVEN `*NODE` followed by `1, 0., 0., 0.`
- WHEN parsed
- THEN node 1 SHALL exist at the origin

### Requirement: Element definition and library
The model SHALL define elements via `*ELEMENT, TYPE=...`, each with an integer element
number and an ordered connectivity list of node numbers. The Phase-1 library SHALL
support the tetrahedral continuum types C3D4 and C3D10; later phases SHALL extend it to
the full reference library (C3D8/C3D20/C3D6/C3D15 and reduced-integration variants,
shell, membrane, beam, truss, plane, axisymmetric, gap, spring, dashpot, and mass
elements). (ref: src/elements.f, src/calinput.f)

#### Scenario: Element with wrong node count
- GIVEN an element card whose connectivity list has the wrong number of nodes for its type
- WHEN parsed
- THEN the program SHALL report an error

#### Scenario: Element references undefined node
- GIVEN an element connectivity referencing a node number that has not been defined
- WHEN the model is assembled
- THEN the program SHALL report an error

### Requirement: 1-D and 2-D element expansion
Beam, shell, and other structural 1-D/2-D elements SHALL be internally expanded into
3-D elements for the solver, governed by their assigned section properties. This is a
later-phase capability and is out of scope for the Phase-1 tetrahedral slice.
(ref: src/gen3delem.f, src/gen3dnor.f)

#### Scenario: Shell expanded to solid
- GIVEN a shell element with a `*SHELL SECTION` thickness
- WHEN the model is prepared for the solver
- THEN the shell SHALL be expanded into a 3-D element using the section thickness

### Requirement: Node and element sets
The model SHALL define named node sets via `*NSET` and element sets via `*ELSET`,
referencing explicit member lists, ranges (`first, last, increment` with the `GENERATE`
parameter), or previously defined sets. Set names SHALL be referenceable by loads,
boundary conditions, sections, and output requests. (ref: src/noset.f, src/elset.f,
src/cascade.c)

#### Scenario: Generated set
- GIVEN `*NSET, NSET=ALL, GENERATE` followed by `1, 100, 1`
- WHEN parsed
- THEN set ALL SHALL contain nodes 1 through 100

### Requirement: Surfaces
The model SHALL define named surfaces via `*SURFACE` as either element-face based
(`TYPE=ELEMENT`) or node based (`TYPE=NODE`), for use in distributed loads, ties, and
contact. Element-face surfaces SHALL resolve to concrete faces of the meshed elements.
(ref: src/surfaces.f)

#### Scenario: Element-face surface
- GIVEN `*SURFACE, NAME=TOP, TYPE=ELEMENT` listing element faces
- WHEN a pressure load references surface TOP
- THEN the load SHALL be applied to the listed element faces

### Requirement: Reference-equivalent model
For a given `test/` deck, the assembled node/element/set/surface model SHALL be
equivalent to what reference CalculiX builds, so downstream solve results match within
the documented numerical tolerance. The model SHALL be inspectable from the Python
bindings (see python-bindings spec).

#### Scenario: Reference-equivalent assembly
- GIVEN a Phase-1 tetrahedral `test/` deck
- WHEN CalculixPP and reference CalculiX build their models
- THEN the node coordinates, element connectivities, and set memberships SHALL agree up to ordering and floating-point tolerance
