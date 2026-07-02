# Crack Propagation

## Purpose

Predicts the growth of cracks under cyclic or static loading by computing stress-intensity
factors along crack fronts and advancing the crack shape over increments, selected by
`*CRACK PROPAGATION`. In CalculixPP the front-geometry math and SIF evaluation are pure
C++20 kernels on **NumPP** containers; crack-shape advance and any front remeshing use the
**CyberCadKernel** geometry kernel (see the `mesh-processing` capability). Field solves feed
through the **ComputeBackend** (CPU default and always available; GPU optional and never
required). The procedure is reachable from the Python bindings. (ref: src/CalculiX.c:1629
crackpropagation, src/crackpropagation.c, src/crackpropagations.f, src/crackrate.f,
src/crackshape.f, src/crackfrd.c)

**Porting Phase:** 5

## Requirements

### Requirement: Crack propagation procedure
A `*CRACK PROPAGATION` step SHALL read a crack definition (crack mesh / front), evaluate
fracture parameters along the crack front from stress results, and advance the crack
incrementally according to a crack-growth law, preserving the reference keyword semantics.
(ref: src/crackpropagation.c)

#### Scenario: Crack advance increment
- GIVEN a `*CRACK PROPAGATION` step over a model with an initial crack
- WHEN an increment is processed
- THEN the solver SHALL compute stress-intensity factors along the front using NumPP kernels and advance the front by the increment dictated by the growth law
- AND the step SHALL be invocable from the Python bindings

### Requirement: Crack length and shape tracking
The procedure SHALL track the evolving crack length and shape, computing crack-front
geometry (local 1-D crack length, smoothing of the front) and delegating shape advance and
front remeshing to CyberCadKernel. (ref: src/cracklength.f, src/cracklength_smoothing.f,
src/crackshape.f)

#### Scenario: Front smoothing
- GIVEN a crack front computed after an advance
- WHEN the next increment is prepared
- THEN the front geometry SHALL be smoothed to avoid irregular jumps before CyberCadKernel updates the crack shape

### Requirement: Crack-growth rate law
The crack advance per cycle SHALL be governed by a configured crack-rate law relating growth
to the stress-intensity range, evaluated by a C++20 kernel. (ref: src/crackrate.f)

#### Scenario: Growth from stress-intensity range
- GIVEN a stress-intensity range computed along the front
- WHEN the crack rate is evaluated
- THEN the per-increment advance SHALL follow the configured crack-growth-rate law

### Requirement: Crack results output
The procedure SHALL write crack-front results (stress-intensity factors and crack length
along the front) to dedicated output for postprocessing. (ref: src/crackfrd.c)

#### Scenario: Write crack results
- GIVEN a completed crack-propagation increment
- WHEN results are written
- THEN the stress-intensity factors and crack length along the front SHALL be output for postprocessing

### Requirement: Reference-result fidelity
Crack-propagation results (stress-intensity factors and crack length per increment) for a given reference `test/` deck SHALL match the reference CalculiX output within the documented numerical tolerance, on the CPU backend with no GPU present. (ref: src/crackfrd.c)

#### Scenario: Match reference crack deck
- GIVEN a reference `*CRACK PROPAGATION` `*.inp` deck run on the CPU ComputeBackend with no GPU toolkit installed
- WHEN the analysis completes
- THEN the computed stress-intensity factors and crack advance SHALL agree with the reference CalculiX results within tolerance
