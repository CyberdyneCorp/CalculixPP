# CFD and Network (1-D Fluid) Analysis

## Purpose

CalculixPP solves fluid behavior in two forms, ported to pure C++20: 3-D
computational fluid dynamics (CFD) on a volume mesh, and 1-D aerodynamic/hydraulic
networks built from fluid-section elements (orifices, pipes, valves, labyrinths)
frequently coupled to thermal analysis. Implicit CFD and network linear solves run
through **NumPP**, with assembly/solve dispatched via the **ComputeBackend** (CPU
always available; GPU optional). Volume meshes and geometry prep reference
**CyberCadKernel** (see mesh-processing). (ref: src/compfluidfem.c,
src/fluidsections.f, src/network*.f, src/radflowload.c, src/flowoutput.f)

**Porting Phase:** 5

## Requirements

### Requirement: 3-D CFD analysis
A `*CFD` step SHALL solve the (incompressible or compressible) Navier-Stokes
equations on a 3-D fluid mesh for velocity, pressure, and (when thermal) temperature,
supporting time-marching to steady state or transient. Field solves SHALL use NumPP
on the ComputeBackend, and the CFD mesh SHALL be prepared via CyberCadKernel
(mesh-processing). The step SHALL be reachable from the Python bindings.

#### Scenario: CFD step
- GIVEN a `*CFD` step with inlet/outlet boundary conditions
- WHEN the step runs
- THEN CalculixPP SHALL output the fluid velocity, pressure, and temperature fields

### Requirement: 1-D fluid networks
CalculixPP SHALL solve gas/liquid networks made of fluid-section elements connected
at nodes, computing mass flow, total pressure, and total temperature throughout the
network. The network solver SHALL handle the nonlinear flow relations of each element
type, with the Newton linear systems solved through NumPP. Fluid-section card
semantics SHALL match the reference.

#### Scenario: Orifice mass flow
- GIVEN a network with an orifice fluid section between two pressure nodes
- WHEN the network is solved
- THEN the mass flow through the orifice SHALL satisfy its characteristic flow equation for the resulting pressure ratio

### Requirement: Coupling to heat transfer
Network analysis SHALL couple to solid heat transfer through forced-convection films,
so gas temperatures and heat exchange with solid walls are solved consistently within
a coupled step (see heat-transfer-analysis).

#### Scenario: Gas temperature drives wall film
- GIVEN a coupled step where a network supplies a forced-convection film on a solid wall
- WHEN the step is solved
- THEN the film temperature SHALL be taken from the network gas temperature and heat SHALL exchange with the wall

### Requirement: Network MPCs
CalculixPP SHALL support `*NETWORK MPC` relations and special network constraints to
model components such as splitters and junctions, enforcing mass-flow conservation
per the constraint.

#### Scenario: Network junction
- GIVEN a `*NETWORK MPC` modeling a flow split at a junction
- WHEN the network is solved
- THEN mass-flow conservation SHALL be enforced at the junction per the constraint

### Requirement: Reference fidelity
CFD and network results SHALL match the reference CalculiX output for the
corresponding `test/` deck within the documented numerical tolerance, on the CPU
ComputeBackend with no GPU present and on any optional GPU backend.

#### Scenario: Regression against reference deck
- GIVEN a reference CFD or network `*.inp` deck with known CalculiX results
- WHEN CalculixPP solves it on the CPU backend
- THEN the fluid fields (or mass flow, pressure, temperature) SHALL agree with the reference within tolerance

### Requirement: Fluid and physical constants (Phase 5)
CalculixPP SHALL parse `*FLUID CONSTANTS` (card `*FLUIDCONSTANTS`, specific heat and dynamic viscosity as functions of temperature), `*SPECIFIC GAS CONSTANT` (card `*SPECIFICGASCONSTANT`, the gas constant R for the ideal-gas law), and `*PHYSICAL CONSTANTS` (card `*PHYSICALCONSTANTS`, absolute-zero offset, Stefan-Boltzmann constant, and gravitational constant), making these constants available to the network/CFD and radiation solves run through NumPP on the ComputeBackend. (ref: src/fluidconstants.f, src/specificgasconstants.f, src/physicalconstants.f)

#### Scenario: Temperature-dependent fluid properties
- GIVEN a `*FLUID CONSTANTS` table of specific heat and viscosity versus temperature
- WHEN the network is solved at a given gas temperature
- THEN the specific heat and viscosity SHALL be interpolated from the table at that temperature

#### Scenario: Ideal-gas relation uses R
- GIVEN a `*SPECIFIC GAS CONSTANT` value R and a `*PHYSICAL CONSTANTS` absolute-zero offset
- WHEN a compressible network element evaluates density from total pressure and temperature
- THEN the ideal-gas relation SHALL use R and the absolute-temperature scale from those cards

### Requirement: Face and flow boundary conditions (Phase 5)
CalculixPP SHALL support face and flow boundary conditions for fluid analysis: `*BOUNDARYF` to prescribe a boundary condition on element faces, `*MASS FLOW` (card `*MASSFLOW`) to prescribe a network mass flow, and `*TRANSFORMF` to define a local coordinate transform on element faces, applied consistently in the network/CFD solves dispatched via the ComputeBackend and reachable from the Python bindings. (ref: src/boundaryfs.f, src/massflows.f, src/transformfs.f)

#### Scenario: Prescribed mass flow in a network
- GIVEN a `*MASS FLOW` card prescribing the mass flow at a network node
- WHEN the network is solved
- THEN the mass flow at that node SHALL be enforced to the prescribed value

#### Scenario: Face boundary condition with local transform
- GIVEN a `*BOUNDARYF` on element faces together with a `*TRANSFORMF` local system on those faces
- WHEN the fluid step is solved
- THEN the face boundary condition SHALL be applied in the local transformed frame
