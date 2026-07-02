# Element Sections and Orientation

## Purpose

Binds materials and geometric/section properties to element sets, and defines the local
material coordinate systems used for anisotropic materials and result output. The
section binding drives which material stiffness (evaluated via NumPP) each element
uses during assembly. (ref: src/solidsections.f, src/shellsections.f, src/beamsections.f,
src/fluidsections.f, src/orientations.f)

**Porting Phase:** 1 — Foundation & vertical slice (`*SOLID SECTION` and `*ORIENTATION` for the linear-static tetrahedral slice; shell, beam, membrane, and fluid sections land in later phases).

## Requirements

### Requirement: Solid sections
The model SHALL assign a material to solid (continuum) elements via
`*SOLID SECTION, ELSET=..., MATERIAL=...`, optionally with an `ORIENTATION`. Every
continuum element used in a solve SHALL resolve to exactly one solid section.
(ref: src/solidsections.f)

#### Scenario: Element with no section
- GIVEN a continuum element belonging to no section
- WHEN the model is assembled
- THEN the program SHALL report an error that the element has no material assigned

### Requirement: Shell and membrane sections
The model SHALL assign properties to shell elements via `*SHELL SECTION` (with a
thickness and optional offset) and to membrane elements via `*MEMBRANE SECTION`; the
thickness SHALL drive the 2-D-to-3-D expansion of the elements. This is a later-phase
capability, out of scope for the Phase-1 tetrahedral slice.
(ref: src/shellsections.f, src/membranesections.f, src/gen3dnor.f)

#### Scenario: Shell thickness
- GIVEN `*SHELL SECTION, ELSET=skin, MATERIAL=al` with a thickness value
- WHEN the model is assembled
- THEN the skin elements SHALL use that thickness for bending and membrane stiffness

### Requirement: Beam sections
The model SHALL assign properties to beam elements via `*BEAM SECTION` (cross-section
shape such as RECT, CIRC, PIPE, BOX with dimensions and section-axis orientation) and
via `*BEAM GENERAL SECTION` (direct stiffness properties). This is a later-phase
capability, out of scope for the Phase-1 tetrahedral slice.
(ref: src/beamsections.f, src/beamgeneralsections.f)

#### Scenario: Rectangular beam
- GIVEN `*BEAM SECTION, SECTION=RECT` with width and height and a section orientation
- WHEN the beam stiffness is built
- THEN the cross-section moments of inertia SHALL be derived from the given dimensions

### Requirement: Fluid sections
The model SHALL assign properties to network/fluid elements via `*FLUID SECTION`,
selecting a fluid-element behavior type (orifice, pipe, valve, etc.) and its geometric
parameters, for use in network (gas/liquid) analysis. This is a later-phase capability
(Phase 5), out of scope for the Phase-1 tetrahedral slice. (ref: src/fluidsections.f)

#### Scenario: Orifice fluid section
- GIVEN `*FLUID SECTION` of an orifice type with its geometric parameters
- WHEN the network is solved
- THEN the element SHALL use the orifice flow characteristic for its mass-flow / pressure relation

### Requirement: Local orientation
The model SHALL define local coordinate systems via `*ORIENTATION` (rectangular or
cylindrical), referenceable from section cards. Anisotropic/orthotropic material
constants and requested local results SHALL be expressed in this system, with the local
axes computed consistently with reference CalculiX. (ref: src/orientations.f)

#### Scenario: Cylindrical orientation
- GIVEN `*ORIENTATION, SYSTEM=CYLINDRICAL` with two points defining the axis
- WHEN material stiffness is evaluated at an integration point
- THEN the local axes SHALL follow the cylindrical system at that point

### Requirement: Discrete and connector element properties (Phase 2)
The model SHALL define the properties of discrete/connector elements via `*MASS` (a concentrated point mass on a node set), `*SPRING` (spring stiffness, optionally nonlinear force-displacement or per-DOF) and `*DASHPOT` (viscous dashpot coefficient). Each card SHALL bind its property to the referenced `ELSET` and contribute the corresponding term (mass, stiffness, or damping) to the assembled matrices during the solve. This is a later-phase capability, out of scope for the Phase-1 tetrahedral slice. (ref: src/masses.f, src/springs.f, src/dashpots.f)

#### Scenario: Point mass contributes to mass matrix
- GIVEN a `*MASS, ELSET=lumped` card with a mass value
- WHEN the mass matrix is assembled for a dynamic or modal procedure
- THEN the given mass SHALL be added to the translational DOFs of the referenced nodes

#### Scenario: Spring stiffness contributes to stiffness matrix
- GIVEN a `*SPRING, ELSET=springs` card with a stiffness value
- WHEN the global stiffness matrix is assembled
- THEN the spring SHALL contribute its stiffness between the connected DOFs

#### Scenario: Dashpot contributes to damping matrix
- GIVEN a `*DASHPOT, ELSET=dampers` card with a viscous coefficient
- WHEN the damping matrix is assembled for a transient procedure
- THEN the dashpot SHALL contribute its coefficient between the connected DOFs

### Requirement: User element and section extensibility (Phase 2)
The model SHALL support user-defined elements via `*USER ELEMENT` (card `*USERELEMENT`) and user-defined sections via `*USER SECTION` (card `*USERSECTION`), exposing a C++20 element/section interface (Python-reachable) that supplies the element stiffness, mass, and internal-force contributions in place of the built-in element routines. The interface SHALL receive the element nodal coordinates, state, and material data, and its contributions SHALL be assembled like those of built-in elements. This is a later-phase capability, out of scope for the Phase-1 tetrahedral slice. (ref: src/userelements.f, src/usersections.f)

#### Scenario: User element supplies its stiffness
- GIVEN a `*USER ELEMENT`/`*USERELEMENT` bound to a C++20 element implementation and a `*USERSECTION` assigning it to an `ELSET`
- WHEN the global stiffness matrix is assembled
- THEN CalculixPP SHALL call the user element to obtain its stiffness contribution and assemble it into the global matrix

#### Scenario: User section reachable from Python
- GIVEN a user-defined section registered through the Python bindings
- WHEN the model referencing that section is assembled
- THEN the Python-reachable interface SHALL be invoked to provide the section's element data

### Requirement: Nodal thickness, normals, and distributions (Phase 2)
The model SHALL support per-node geometric data and spatially-varying element data via `*NODAL THICKNESS` (card `*NODALTHICKNESS`) assigning thickness to individual nodes of shells/beams, `*NORMAL` defining user-specified nodal normals for shells/beams, and `*DISTRIBUTION` with its governing `*DISTRIBUTION TABLE` (card `*DISTRIBUTIONTABLE`) supplying spatially-varying per-element data such as element orientation. The `*DISTRIBUTIONTABLE` SHALL define the type and units of the distributed quantity, and nodal thickness/normals SHALL override the section defaults during the 2-D-to-3-D expansion. This is a later-phase capability, out of scope for the Phase-1 tetrahedral slice. (ref: src/nodalthicknesses.f, src/normals.f, src/distributions.f, src/distributiontables.f, src/gen3dnor.f)

#### Scenario: Nodal thickness overrides section thickness
- GIVEN a shell with `*NODALTHICKNESS` values at its nodes
- WHEN the elements are expanded to 3-D
- THEN the per-node thickness SHALL be used instead of the uniform section thickness

#### Scenario: User-defined nodal normal
- GIVEN a `*NORMAL` card specifying a normal direction at a shell/beam node
- WHEN the element local frame is built
- THEN the given normal SHALL be used instead of the geometrically computed normal

#### Scenario: Per-element orientation from a distribution
- GIVEN a `*DISTRIBUTION` of element orientations governed by a `*DISTRIBUTIONTABLE` defining its type and units
- WHEN material stiffness is evaluated at an integration point of an element in the distribution
- THEN the element's local axes SHALL follow its per-element entry in the distribution
