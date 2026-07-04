# High-Cycle Fatigue Specification

## MODIFIED Requirements

### Requirement: Request HCF evaluation over prior results
A `*HCF` step SHALL request a stress-life high-cycle-fatigue evaluation over the stress field of the preceding stress-producing model, selecting the uniaxial-equivalent stress criterion (signed or plain von Mises) that reduces each point's stress tensor to a scalar amplitude.

The shipped physics is a stress-life (Basquin S-N) evaluation: the recovered stress field is treated as the cyclic stress amplitude at each point. `CRITERION=` selects the scalar reduction (`SIGNED-VON-MISES` default, or `VON-MISES`); any other criterion is rejected. This is NOT the CalculiX crack-growth cumulative HCF path (`hcfs.f`/`combilcfhcf.f`), which is a documented follow-on. (ref: src/hcfs.f as the contrasting crack-growth path)

#### Scenario: HCF card triggers evaluation
- GIVEN a stress-producing model (an elastic material and loads) and an `*HCF` card
- WHEN the analysis reaches the HCF request
- THEN the criterion SHALL be evaluated over the recovered stress field on the CPU ComputeBackend
- AND the evaluation SHALL be invocable from the Python bindings

#### Scenario: Criterion selection
- GIVEN an `*HCF, CRITERION=VON-MISES` card
- WHEN the deck is parsed
- THEN the plain von Mises amplitude SHALL be used at each point
- AND an unrecognized `CRITERION=` value SHALL raise a parse error

### Requirement: Fatigue criterion produces worst-case location and life
The `*HCF` evaluation SHALL invert the material S-N (Basquin) curve `S_a = a·N^b` at each result location to a cycles-to-failure `N = (S_a/a)^(1/b)`, and report the critical (worst-case) location — the point of largest stress amplitude, hence smallest life — together with its estimated fatigue life.

The per-node scalar amplitude `S_a` comes from the selected criterion; the material's `*FATIGUE` card supplies the Basquin coefficient `a` and exponent `b` (`b < 0`). The report carries the per-node life and amplitude plus the worst-case node id, location, amplitude, and life. (ref: Basquin stress-life S_a = a·N^b)

#### Scenario: Critical location reported
- GIVEN an `*HCF` evaluation over a stress field with a well-defined maximum-amplitude point
- WHEN the criterion has been applied across all result locations
- THEN the location with the largest stress amplitude (smallest life) SHALL be identified as the critical location
- AND its estimated life `N = (S_a/a)^(1/b)` SHALL be reported via the results output

### Requirement: Require a valid preceding results source
The `*HCF` request SHALL require a valid stress-producing model and an S-N material curve, and SHALL raise a diagnostic error naming the missing source rather than produce output when either is absent.

#### Scenario: Missing S-N curve is an error
- GIVEN an `*HCF` deck whose materials carry no `*FATIGUE` (S-N) curve
- WHEN the HCF evaluation is attempted
- THEN the solver SHALL raise a clear diagnostic error identifying the missing S-N curve
- AND SHALL NOT emit fatigue results

#### Scenario: Missing stress source is an error
- GIVEN an `*HCF` deck with no elastic material to recover a stress field from
- WHEN the HCF evaluation is attempted
- THEN the solver SHALL raise a clear diagnostic error identifying the missing stress source
- AND SHALL NOT emit fatigue results

### Requirement: Reference-result fidelity
The `*HCF` life estimate SHALL match the closed-form Basquin inversion `N = (S_a/a)^(1/b)` for a point of known stress amplitude within the documented numerical tolerance, on the CPU backend with no GPU present.

This replaces the CalculiX crack-growth reference-fidelity requirement for the HCF procedure: the CalculiX `*HCF` output is a crack-growth cumulative-damage quantity (different physics) and would not match a stress-life estimate, so fidelity is asserted against the closed-form Basquin curve instead. (ref: Basquin S_a = a·N^b; src/hcfs.f as the contrasting crack-growth path)

#### Scenario: Match the closed-form Basquin life
- GIVEN a uniaxial specimen at a known stress amplitude `S` and an S-N curve `(a, b)` run on the CPU ComputeBackend with no GPU toolkit installed
- WHEN the fatigue evaluation completes
- THEN the reported life at that point SHALL equal `(S/a)^(1/b)` within `1e-6` relative tolerance
- AND the worst-case location SHALL be the highest-amplitude point
