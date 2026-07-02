# Multi-Point Constraints and Connections

## Purpose

Defines relationships that tie degrees of freedom of different nodes together: linear
equations, multi-point constraints (MPCs), rigid bodies, couplings, ties between
surfaces, and cyclic symmetry. In CalculixPP these constraints are eliminated or
enforced during assembly; the dependent-DOF elimination linear solves are performed via
NumPP (replacing the Fortran cascade). Results SHALL match reference CalculiX for the
corresponding test deck within documented tolerance. (ref: src/equations.f,
src/rigidbodys.f, src/couplings.f, src/ties.f, src/cyclicsymmetrymodels.f, src/cascade.c)

**Porting Phase:** 2 — Nonlinear statics & materials (equations, MPCs, rigid bodies,
couplings, ties); cyclic symmetry lands in Phase 4 — Dynamics & eigenproblems (noted per
requirement).

## Requirements

### Requirement: Linear equations
The model SHALL support `*EQUATION`, expressing a linear homogeneous relation among
selected nodal degrees of freedom. The first term's DOF SHALL be the dependent
(eliminated) term, eliminated during assembly using NumPP linear algebra.
(ref: src/equations.f)

#### Scenario: Dependent DOF also constrained by SPC
- GIVEN a DOF that is the dependent term of an `*EQUATION`
- AND the same DOF is also given a `*BOUNDARY` condition (see loads-and-boundary-conditions)
- WHEN the model is assembled
- THEN CalculixPP SHALL report an over-constraint error

### Requirement: Multi-point constraints
The model SHALL support analytic MPC types via `*MPC` (e.g. PLANE, STRAIGHT, BEAM) and
user-defined MPCs, expanded internally into linear equations. (ref: src/calinput.f,
src/cascade.c)

#### Scenario: PLANE MPC
- GIVEN a `*MPC, PLANE` applied to a set of nodes
- WHEN the model is assembled
- THEN the listed nodes SHALL be constrained to remain in a plane defined by the constraint

### Requirement: Rigid bodies
The model SHALL support `*RIGID BODY` tying a node set (or element set) to a reference
node so that the set moves as a rigid body. (ref: src/rigidbodys.f)

#### Scenario: Rigid body motion
- GIVEN a `*RIGID BODY` tying a node set to a reference node
- WHEN the reference node is displaced or rotated
- THEN every node in the set SHALL move rigidly with the reference node
- AND the response SHALL match reference CalculiX within tolerance

### Requirement: Distributing and kinematic couplings
The model SHALL support `*COUPLING` with `*KINEMATIC` (constrains selected DOFs of a
surface to a reference node rigidly) and `*DISTRIBUTING` / `*DISTRIBUTING COUPLING`
(distributes loads/motion of a reference node over a surface by weights).
(ref: src/couplings.f, src/distributingcouplings.f)

#### Scenario: Load through a distributing coupling
- GIVEN a `*COUPLING, *DISTRIBUTING` from a reference node to a surface
- WHEN a load is applied at the reference node (see loads-and-boundary-conditions)
- THEN the load SHALL be distributed across the coupled surface nodes by their weights

### Requirement: Surface ties
The model SHALL support `*TIE` to bond two surfaces, including mesh-independent (mortar)
tied contact and cyclic-symmetry ties. A tolerance parameter SHALL control
node-to-surface association. (ref: src/ties.f, src/tiefaccont.f)

#### Scenario: Tied dissimilar meshes
- GIVEN a `*TIE` between two non-matching meshed surfaces
- WHEN the model is assembled
- THEN displacements SHALL be made continuous across the interface without requiring matching nodes

### Requirement: Cyclic symmetry (Phase 4)
The model SHALL support `*CYCLIC SYMMETRY MODEL` defining a sector with master and slave
surfaces and a number of sectors, enabling cyclic-symmetric static and modal analyses,
with `*SELECT CYCLIC SYMMETRY MODES` selecting nodal diameters. The cyclic-symmetric
eigenproblem SHALL be solved via the NumPP eigensolver (replacing ARPACK).
(ref: src/cyclicsymmetrymodels.f, src/arpackcs.c)

#### Scenario: Modal analysis of a cyclic sector
- GIVEN a `*CYCLIC SYMMETRY MODEL` and a `*FREQUENCY` step
- WHEN the eigenmodes are computed
- THEN CalculixPP SHALL compute eigenfrequencies per selected nodal diameter using the cyclic-symmetric eigensolver
