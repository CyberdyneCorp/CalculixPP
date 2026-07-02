# High-Cycle Fatigue

## Purpose

Estimates high-cycle-fatigue (HCF) life and critical locations from an already-computed
stress result — typically a cyclic, modal, or steady-state dynamic solution — by applying
a fatigue/damage criterion over the stored field. In CalculixPP the criterion evaluation
runs as C++20 kernels over **NumPP** containers dispatched through the **ComputeBackend**
(CPU default and always available; CUDA/OpenCL/Metal optional and never required), and the
evaluation is reachable from the Python bindings. It consumes results produced under
**Results and Output** from a preceding **Eigenvalue Analysis (Frequency, Buckling, Complex
Frequency)** or **Dynamic Analysis** run. (ref: src/CalculiX.c HCF path, src/hcf.f)

**Porting Phase:** 5 — Advanced physics

## Requirements

### Requirement: Request HCF evaluation over prior results
The solver SHALL support `*HCF` to request a high-cycle-fatigue evaluation over previously computed stress results, selecting the results source and the fatigue/damage criterion to apply.
(ref: src/hcf.f)

#### Scenario: HCF card triggers evaluation
- GIVEN a completed stress result and an `*HCF` card naming a fatigue criterion
- WHEN the analysis reaches the HCF request
- THEN the criterion SHALL be evaluated over the stored stress field through the ComputeBackend
- AND the evaluation SHALL be invocable from the Python bindings

### Requirement: Fatigue criterion produces worst-case location and life
The `*HCF` evaluation SHALL apply the selected fatigue/damage criterion at each result location and report the critical (worst-case) location together with its estimated fatigue life or safety measure.
(ref: src/hcf.f)

#### Scenario: Critical location reported
- GIVEN an `*HCF` evaluation over a cyclic or steady-state stress result
- WHEN the criterion has been applied across all result locations
- THEN the location with the most severe fatigue damage SHALL be identified
- AND its estimated life or safety factor SHALL be reported via the results output

### Requirement: Require a valid preceding results source
The `*HCF` request SHALL require a valid, completed preceding results source and SHALL raise a diagnostic error rather than produce output when no such source is available.
(ref: src/hcf.f)

#### Scenario: Missing results source is an error
- GIVEN an `*HCF` card whose named results source has not been computed or is absent
- WHEN the HCF evaluation is attempted
- THEN the solver SHALL raise a clear diagnostic error identifying the missing source
- AND SHALL NOT emit fatigue results

### Requirement: Reference-result fidelity
HCF outputs (critical location and estimated life/safety measure) for a reference `test/` deck SHALL match the reference CalculiX output within the documented numerical tolerance, on the CPU backend with no GPU present.
(ref: src/CalculiX.c HCF path, src/hcf.f)

#### Scenario: Match reference HCF deck
- GIVEN a reference `*HCF` `*.inp` deck run on the CPU ComputeBackend with no GPU toolkit installed
- WHEN the fatigue evaluation completes
- THEN the reported critical location and life estimate SHALL agree with the reference CalculiX results within tolerance
