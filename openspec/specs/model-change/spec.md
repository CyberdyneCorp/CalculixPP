# Model Change

## Purpose

Activates or deactivates parts of the model between steps — element birth/death and
contact-pair activation — so a single analysis can represent excavation, additive
manufacturing layer deposition, or staged assembly. In CalculixPP the model-change
bookkeeping (which elements and contact pairs contribute to the current step) drives the
assembly performed by the **ComputeBackend** over **NumPP** containers (CPU default and
always available; CUDA/OpenCL/Metal optional and never required), and every operation is
reachable from the Python bindings. See **Contact** and **Static Analysis** for the step
machinery it modifies. (ref: src/CalculiX.c model-change path, src/modelchanges.f)

**Porting Phase:** 3 — Thermal & coupled / general nonlinear

## Requirements

### Requirement: Deactivate an element set
The solver SHALL support `*MODEL CHANGE` (card `*MODELCHANGE`) with `TYPE=ELEMENT, REMOVE` to deactivate an element set for the current and subsequent steps so those elements contribute no stiffness, mass, or stress to the assembled system.
(ref: src/modelchanges.f)

#### Scenario: Removed elements carry no load
- GIVEN a converged model and a `*MODEL CHANGE, TYPE=ELEMENT, REMOVE` referencing an element set in a new step
- WHEN the step is assembled through the ComputeBackend
- THEN the removed elements SHALL contribute no stiffness or stress to the NumPP system
- AND the operation SHALL be invocable from the Python bindings

### Requirement: Reactivate an element set
The solver SHALL support `*MODEL CHANGE` with `TYPE=ELEMENT, ADD` to reactivate a previously removed element set, reintroducing it strain-free so its stress state begins from the deformed configuration at reactivation.
(ref: src/modelchanges.f)

#### Scenario: Strain-free reactivation
- GIVEN an element set previously removed with `*MODEL CHANGE, TYPE=ELEMENT, REMOVE`
- WHEN a later step issues `*MODEL CHANGE, TYPE=ELEMENT, ADD` for that set
- THEN the reactivated elements SHALL rejoin the assembly with zero initial strain relative to the current deformed geometry
- AND SHALL begin contributing stiffness and stress from that step onward

### Requirement: Activate or deactivate a contact pair
The solver SHALL support `*MODEL CHANGE` with `TYPE=CONTACT PAIR` (`ADD` or `REMOVE`) to enable or disable a named contact pair between steps, so the contact search and constraint assembly described in **Contact** run only while the pair is active.
(ref: src/modelchanges.f)

#### Scenario: Contact pair disabled between steps
- GIVEN an active contact pair and a `*MODEL CHANGE, TYPE=CONTACT PAIR, REMOVE` for it in a new step
- WHEN the step is solved
- THEN no contact search, force, or constraint SHALL be assembled for that pair through the ComputeBackend
- AND a later `*MODEL CHANGE, TYPE=CONTACT PAIR, ADD` SHALL re-enable it

### Requirement: Reference-result fidelity
Results for a reference `test/` deck using `*MODEL CHANGE` SHALL match the reference CalculiX output within the documented numerical tolerance, on the CPU backend with no GPU present.
(ref: src/CalculiX.c model-change path, src/modelchanges.f)

#### Scenario: Match reference staged deck
- GIVEN a reference excavation or staged-assembly `*.inp` deck run on the CPU ComputeBackend with no GPU toolkit installed
- WHEN the analysis completes across all steps
- THEN the displacements, stresses, and reaction forces SHALL agree with the reference CalculiX results within tolerance
