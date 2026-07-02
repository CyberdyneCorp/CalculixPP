## ADDED Requirements

### Requirement: Spatial proximity search
The engine SHALL find candidate contacting slave-node/master-face pairs each increment using a spatial acceleration structure, bounded by a configurable search/adjust distance, so contact detection scales to large surfaces without an all-pairs test. (ref: src/gencontelem_n2f.f, src/gencontelem_f2f.f)

#### Scenario: Candidate pairs within search distance
- GIVEN slave and master contact surfaces and a search distance
- WHEN proximity search runs for an increment
- THEN only master faces within the search distance of each slave entity SHALL be returned as candidates

#### Scenario: Large surfaces avoid all-pairs cost
- GIVEN contact surfaces with many faces
- WHEN search runs
- THEN the engine SHALL use the spatial structure rather than an O(n²) comparison of every slave against every master face

### Requirement: Node-to-surface projection
The engine SHALL project each slave node onto its nearest candidate master face, computing the signed gap and a local normal/tangent basis at the projection point. (ref: src/gencontelem_n2f.f)

#### Scenario: Slave node projected to master face
- GIVEN a slave node and its candidate master faces
- WHEN projection runs
- THEN the engine SHALL return the closest projection point, the signed normal gap, and the local normal/tangent frame

### Requirement: Surface-to-surface mortar integration
The engine SHALL, for the surface-to-surface formulation, integrate contact contributions over the overlapping segments of master and slave faces using a dual or standard Lagrange basis, so non-matching meshes transfer pressure without over-stiffening. (ref: src/contactmortar.c, src/gencontelem_f2f.f)

#### Scenario: Non-matching mesh mortar coupling
- GIVEN slave and master surfaces with non-matching discretizations
- WHEN mortar integration runs
- THEN contact contributions SHALL be integrated over the overlap segments using the selected (dual/standard) basis

### Requirement: Gap and pressure-overclosure evaluation
The engine SHALL evaluate the normal gap/penetration and derive the normal contact pressure from the active pressure-overclosure law supplied by the contact surface behavior, and SHALL evaluate the tangential stick/slip state from the friction model. (ref: src/surfacebehaviors.f, src/frictions.f)

#### Scenario: Normal pressure from overclosure law
- GIVEN a detected penetration and a hard/exponential/linear pressure-overclosure law
- WHEN the contact state is evaluated
- THEN the normal pressure SHALL follow that law for the current gap

#### Scenario: Tangential stick to slip transition
- GIVEN a contacting point under increasing tangential load with Coulomb friction
- WHEN the tangential traction reaches the friction limit
- THEN the point SHALL transition from stick to slip

### Requirement: Contact operator generation through the ComputeBackend
The engine SHALL generate the contact stiffness/residual contributions and assemble them into the global system through the ComputeBackend / NumPP spine, consistent with the nonlinear-solution-control iteration, so contact is just another contribution to the tangent and residual. (ref: src/contactmortar.c)

#### Scenario: Contact contributes to tangent and residual
- GIVEN an active contact set within a Newton iteration
- WHEN the increment is assembled
- THEN the contact contributions SHALL be added to the global tangent and residual via the ComputeBackend and solved through NumPP

### Requirement: Shared search for thermal contact
The engine SHALL expose the current contact state (pressure, clearance, contact area) to the thermal model so gap conductance and gap heat generation reuse the same search, giving coupled thermomechanical contact a single consistent interface geometry. (ref: src/radflowload.c, heat-transfer-analysis)

#### Scenario: Gap state drives thermal conductance
- GIVEN a coupled thermomechanical contact step
- WHEN the interface is evaluated
- THEN the same gap/pressure state SHALL feed the thermal gap conductance and gap heat generation (see contact, heat-transfer-analysis) without a second search
