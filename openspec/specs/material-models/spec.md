# Material Models

## Purpose

Defines the constitutive behavior of materials assigned to elements: elasticity,
plasticity, thermal/physical properties, density, and specialized models. A material is
introduced by `*MATERIAL, NAME=...` and described by the property cards that follow.
CalculixPP integrates constitutive laws in C++20 kernels (replacing the Fortran element
routines) and exposes a C++ material-model interface (Python-reachable) in place of the
Fortran UMAT. Results SHALL match reference CalculiX for the corresponding test deck
within documented tolerance. (ref: src/materials.f, src/elastics.f, src/plastics.f,
src/umat.f)

**Porting Phase:** 1 — Foundation & vertical slice (linear elasticity + density); plasticity,
hyperelasticity, creep, and user materials land in Phase 2 — Nonlinear statics & materials
(noted per requirement).

## Requirements

### Requirement: Named materials and section assignment
The model SHALL define materials via `*MATERIAL, NAME=...`. A material SHALL be
referenced from a section card (e.g. `*SOLID SECTION, MATERIAL=...`) to bind it to a
set of elements, using the same Abaqus semantics as reference CalculiX.
(ref: src/materials.f)

#### Scenario: Section references undefined material
- GIVEN a section card referencing a material name that was never defined
- WHEN the model is assembled
- THEN CalculixPP SHALL report an error identifying the missing material

### Requirement: Linear elasticity (Phase 1)
A material SHALL define linear elasticity via `*ELASTIC`, supporting `ISO` (isotropic,
Phase 1), and `ORTHO` (orthotropic) / `ANISO` (fully anisotropic) types, with
temperature-dependent constants. The isotropic constitutive matrix SHALL be evaluated
by a C++20 kernel on the CPU ComputeBackend (always available, no GPU required).
(ref: src/elastics.f, src/anisotropic.f)

#### Scenario: Isotropic elastic stress-strain (Phase 1)
- GIVEN an `*ELASTIC, TYPE=ISO` material with Young's modulus and Poisson's ratio
- WHEN a linear `*STATIC` step is solved (see loads-and-boundary-conditions)
- THEN the computed displacements SHALL match reference CalculiX for the deck within tolerance

#### Scenario: Temperature-dependent elasticity
- GIVEN `*ELASTIC` with multiple lines each ending in a temperature
- WHEN a material point is evaluated at a temperature between table entries
- THEN the elastic constants SHALL be linearly interpolated

### Requirement: Density and physical properties
A material SHALL define `*DENSITY` (Phase 1), and `*EXPANSION` (thermal expansion),
`*CONDUCTIVITY`, and `*SPECIFIC HEAT`; these properties SHALL be consumed by the
relevant procedures (dynamics/mass needs density; heat transfer needs conductivity and
specific heat). Thermal expansion couples to `*TEMPERATURE` loading (see
loads-and-boundary-conditions). (ref: src/densities.f, src/expansions.f,
src/conductivities.f, src/specificheats.f)

#### Scenario: Density parsed and stored (Phase 1)
- GIVEN a `*MATERIAL` with a `*DENSITY` card
- WHEN the model is assembled
- THEN the density SHALL be stored on the material and available to mass/body-load computation

#### Scenario: Dynamic step without density
- GIVEN a dynamic or modal procedure and a material lacking a `*DENSITY` card
- WHEN the mass matrix is assembled
- THEN CalculixPP SHALL report an error

### Requirement: Plasticity and hardening (Phase 2)
A material SHALL define rate-independent plasticity via `*PLASTIC` (isotropic and
kinematic hardening) with a yield-stress / plastic-strain table, and cyclic hardening
via `*CYCLIC HARDENING`. Plastic response SHALL be integrated at each integration point
by a C++20 constitutive kernel during nonlinear analysis. (ref: src/plastics.f,
src/nonlingeo.c)

#### Scenario: Yielding beyond the elastic limit
- GIVEN a `*PLASTIC` material with a yield-stress / plastic-strain table
- WHEN the stress at an integration point exceeds the current yield stress
- THEN plastic strain SHALL accumulate and the stress SHALL follow the hardening curve
- AND results SHALL match reference CalculiX for the deck within tolerance

### Requirement: Extended constitutive models and user materials (Phase 2)
The material library SHALL additionally support `*HYPERELASTIC`, `*HYPERFOAM`, `*CREEP`,
`*VISCO` time integration, `*DEFORMATION PLASTICITY`, `*MOHR COULOMB`, and
`*DAMAGE INITIATION` (damage onset criteria), and SHALL expose a C++20 material-model
interface for user-defined behavior via `*USER MATERIAL` (replacing the Fortran UMAT).
The interface SHALL be reachable from the Python bindings.
(ref: src/hyperelastics.f, src/creeps.f, src/usermaterials.f, src/umat.f)

#### Scenario: User-defined material via C++ interface
- GIVEN a `*USER MATERIAL` bound to a C++20 material-model implementation
- WHEN an integration point is evaluated
- THEN CalculixPP SHALL call the user model to obtain the stress and the tangent stiffness

#### Scenario: Large-strain hyperelastic response
- GIVEN a `*HYPERELASTIC` material under large deformation
- WHEN the element response is evaluated
- THEN stress SHALL be derived from the strain-energy potential rather than a linear modulus

#### Scenario: Damage initiation criterion
- GIVEN a material with a `*DAMAGE INITIATION` criterion
- WHEN an integration point reaches the criterion during the solution
- THEN CalculixPP SHALL flag damage onset at that point for output and any coupled degradation

### Requirement: User-material state and rate dependence (Phase 2)
A material SHALL declare the number of solution-dependent internal state variables via `*DEPVAR` for a `*USER MATERIAL`, and MAY define strain-rate-dependent yield scaling via `*RATEDEPENDENT` (card `*RATEDEPENDENT`) for plasticity. The state-variable storage SHALL be allocated per integration point and passed to the C++20 material-model interface (Python-reachable) on each evaluation, and the rate-dependent factor SHALL scale the yield stress as a function of equivalent plastic strain rate. (ref: src/depvars.f, src/ratedependents.f, src/umat.f)

#### Scenario: State variables allocated for user material
- GIVEN a `*USER MATERIAL` with a `*DEPVAR` card declaring N internal state variables
- WHEN an integration point bound to that material is evaluated
- THEN CalculixPP SHALL provide N per-point state variables to the C++20 material model and persist their updated values between increments

#### Scenario: Rate-dependent yield scaling
- GIVEN a `*PLASTIC` material with a `*RATEDEPENDENT` card giving a yield ratio versus equivalent plastic strain rate
- WHEN the plastic response is integrated under a nonzero strain rate
- THEN the yield stress SHALL be scaled by the interpolated rate-dependent factor
- AND results SHALL match reference CalculiX for the deck within tolerance

### Requirement: Additional constitutive data cards (Phase 2)
A material SHALL support supplementary constitutive data cards: `*MOHR COULOMB HARDENING` (card `*MOHRCOULOMBHARDENING`) providing the hardening curve for a `*MOHR COULOMB` material, `*INITIAL STRAIN INCREASE` (card `*INITIALSTRAININCREASE`) applying an increment to the initial (eigen)strain, and `*VALUES AT INFINITY` (card `*VALUESATINFINITY`) giving the long-term (t→∞) elastic constants used by `*VISCO` viscoelasticity. Each card SHALL be parsed only in conjunction with its governing material card and consumed by the corresponding C++20 constitutive kernel. (ref: src/mohrcoulombhardenings.f, src/initialstrainincreases.f, src/valuesatinfinity.f)

#### Scenario: Mohr-Coulomb hardening curve
- GIVEN a `*MOHR COULOMB` material with a `*MOHRCOULOMBHARDENING` hardening table
- WHEN the material yields during a nonlinear step
- THEN the cohesion SHALL evolve along the given hardening curve

#### Scenario: Initial strain increment applied
- GIVEN a material with an `*INITIALSTRAININCREASE` card
- WHEN the initial (eigen)strain field is established for the step
- THEN the specified increment SHALL be added to the initial strain used in the stress computation

#### Scenario: Long-term elastic constants for viscoelasticity
- GIVEN a `*VISCO` material with a `*VALUESATINFINITY` card of long-term elastic constants
- WHEN the time-dependent relaxation is integrated toward t→∞
- THEN the response SHALL approach the stiffness defined by the values-at-infinity constants
