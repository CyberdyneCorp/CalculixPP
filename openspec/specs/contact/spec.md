# Contact

## Purpose

Models interaction between surfaces that may open and close during a nonlinear step,
transferring normal pressure and tangential friction using node-to-surface or
surface-to-surface (mortar) formulations. In CalculixPP the contact search, mortar
integration, and constraint assembly are pure C++20 kernels built on **NumPP** dense/sparse
containers, assembled and solved through the **ComputeBackend** abstraction (CPU default
and always available; CUDA/OpenCL/Metal optional and never required). All contact operations
are reachable from the Python bindings. (ref: src/contactpairs.f, src/contactmortar.c,
src/gencontelem_n2f.f, src/gencontelem_f2f.f, src/surfacebehaviors.f, src/frictions.f)

**Porting Phase:** 3

## Requirements

### Requirement: Contact pairs
The solver SHALL define contact between a slave and master surface via `*CONTACT PAIR`
referencing a `*SURFACE INTERACTION`, selecting the formulation via `TYPE=NODE TO SURFACE`
or `TYPE=SURFACE TO SURFACE` (mortar), preserving the reference keyword semantics.
(ref: src/contactpairs.f)

#### Scenario: Surfaces in contact
- GIVEN a `*CONTACT PAIR` between a slave and master surface
- WHEN the nonlinear step detects penetration
- THEN a normal contact force SHALL be assembled through the ComputeBackend to resist penetration according to the pressure-overclosure relation
- AND the operation SHALL be invocable from the Python bindings

### Requirement: Surface behavior (normal)
A `*SURFACE INTERACTION` SHALL define the normal pressure-overclosure relation via
`*SURFACE BEHAVIOR` (hard, exponential, or linear/tied), evaluated by C++20 kernels.
(ref: src/surfacebehaviors.f)

#### Scenario: Hard pressure-overclosure
- GIVEN a `*SURFACE BEHAVIOR` of hard type
- WHEN slave and master surfaces touch
- THEN pressure SHALL build up to resist penetration with (near-)zero allowed overclosure

### Requirement: Friction
Tangential behavior SHALL be defined via `*FRICTION` with a Coulomb coefficient and a
stick stiffness, transferring shear up to the friction limit, with the stick/slip state
resolved per contact point. (ref: src/frictions.f)

#### Scenario: Stick to slip transition
- GIVEN a contact pair with `*FRICTION`
- WHEN the tangential traction reaches the friction coefficient times the normal pressure
- THEN the contact point SHALL transition from stick to slip

### Requirement: Mortar (surface-to-surface) contact
Surface-to-surface contact SHALL use a mortar formulation that integrates contact
constraints over slave faces, supporting dual or standard Lagrange shape functions for
accurate pressure on non-matching meshes. The mortar coupling matrices SHALL be built with
NumPP sparse containers and assembled via the ComputeBackend. (ref: src/contactmortar.c,
src/multimortar.c, src/mortar.h)

#### Scenario: Pressure on non-matching mesh
- GIVEN a `TYPE=SURFACE TO SURFACE` contact pair between non-matching meshes
- WHEN the step converges
- THEN the contact pressure SHALL be integrated over slave faces using the mortar formulation

### Requirement: Contact modifiers and output
The solver SHALL support `*CLEARANCE`, `*CONTACT DAMPING`, `*CHANGE FRICTION`,
`*CHANGE SURFACE BEHAVIOR`, and `*CHANGE CONTACT TYPE` to modify contact between steps,
and SHALL emit contact results (stress, slip, opening, forces) via `*CONTACT FILE`,
`*CONTACT PRINT`, and `*CONTACT OUTPUT`. (ref: src/calinput.f, src/clearances.f,
src/mortar_postfrd.c)

#### Scenario: Adjust initial clearance and write results
- GIVEN a `*CLEARANCE` setting a uniform initial gap and a `*CONTACT FILE` request for contact stress (CSTR)
- WHEN the step starts and an output increment completes
- THEN the contact gap SHALL be initialized to the specified value regardless of mesh geometry
- AND the contact stresses SHALL be written for postprocessing

### Requirement: Reference-result fidelity
Contact results (contact stress, slip, opening, and reaction forces) for a given reference `test/` deck SHALL match the reference CalculiX output within the documented numerical tolerance, on the CPU backend with no GPU present. (ref: src/mortar_postfrd.c)

#### Scenario: Match reference contact deck
- GIVEN a reference contact `*.inp` deck run on the CPU ComputeBackend with no GPU toolkit installed
- WHEN the analysis completes
- THEN the contact stresses and forces SHALL agree with the reference CalculiX results within tolerance

### Requirement: Thermal contact conductance and gap heat generation (Phase 3)
A `*SURFACE INTERACTION` SHALL support thermal contact via `*GAP CONDUCTANCE` (card `*GAPCONDUCTANCE`) defining heat conductance across the contact gap as a function of contact pressure and/or clearance, and `*GAP HEAT GENERATION` (card `*GAPHEATGENERATION`) defining heat generated in the gap (e.g. frictional or electrical), with the resulting gap heat flux assembled through the ComputeBackend and reachable from the Python bindings. (ref: src/gapconductances.f, src/gapheatgenerations.f)

#### Scenario: Conductance varies with contact pressure
- GIVEN a `*GAP CONDUCTANCE` table giving conductance versus contact pressure
- WHEN two surfaces in a coupled thermal-mechanical step touch under pressure
- THEN the heat conductance across the gap SHALL be interpolated from the table at the current contact pressure
- AND the corresponding gap heat flux SHALL be assembled through the ComputeBackend

#### Scenario: Frictional heat in the gap
- GIVEN a `*GAP HEAT GENERATION` card on a contact pair
- WHEN the surfaces slide under friction during the step
- THEN heat SHALL be generated in the contact gap per the card and distributed to the contacting surfaces
