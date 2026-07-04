# Eigenvalue Analysis (Frequency, Buckling, Complex Frequency) Specification

## MODIFIED Requirements

### Requirement: Complex frequency analysis
A `*COMPLEX FREQUENCY` step SHALL compute complex eigenmodes built from the real eigenmodes of a preceding `*FREQUENCY` step, accounting for proportional (Rayleigh / modal) damping, and report each mode's damped natural frequency, damping ratio, decay rate, and stability.

The shipped physics is proportional-damping subspace-projected complex modes (ABAQUS-style): the real damping operator `C = αM + βK` (and/or `*MODAL DAMPING` ratios) is reduced onto the real modal basis to the quadratic problem `(λ²I + λC_r + Λ)q = 0`, linearized to a real `2·nev` companion pencil, and solved with the dense NumPP `eig()`. This is NOT the CalculiX `*COMPLEX FREQUENCY, CORIOLIS` gyroscopic operator (skew reduced matrix, `i·ω·G_r` coupling, rotor whirl / Campbell split, driven by a rotation body load); the true `CORIOLIS` operator and the `FLUTTER` complex applied-force path are documented follow-ons. (ref: NumPP linalg.hpp `eig`; CalculiX complexfreq.c / coriolissolve.f / complexfrequencys.f as the contrasting gyroscopic path)

#### Scenario: Damped complex modes from a proportional-damping step
- GIVEN a `*FREQUENCY` step followed by `*COMPLEX FREQUENCY` with a proportional `*DAMPING` (Rayleigh) or `*MODAL DAMPING` card
- WHEN the step runs
- THEN the solver SHALL output complex eigenvalues whose imaginary part gives the damped frequency and whose real part gives the decay rate, plus a per-mode damping ratio and a stability flag (`stable` when `ζ > 0`)

#### Scenario: Missing eigenmode source
- GIVEN a `*COMPLEX FREQUENCY` step with no preceding `*FREQUENCY` step in the job
- WHEN the job runs
- THEN the program SHALL report an error that the required eigenmodes are unavailable

### Requirement: Reference-result fidelity
For every real eigenvalue `test/` deck (frequency, buckling), CalculixPP eigenfrequencies and buckling factors SHALL match reference CalculiX output within a documented relative tolerance; for the proportional-damping `*COMPLEX FREQUENCY` procedure, complex eigenvalues (damped frequency and damping ratio) SHALL instead match the analytical closed form within a documented relative tolerance, because CalculiX `*COMPLEX FREQUENCY` solves the different (gyroscopic) eigenproblem and its output would not match. All fidelity checks are validated through the pytest bindings.

#### Scenario: Eigenfrequencies match reference
- GIVEN a reference CalculiX `*FREQUENCY` test deck and its expected eigenfrequencies
- WHEN CalculixPP runs the same deck
- THEN each computed eigenfrequency SHALL match the reference within the documented relative tolerance

#### Scenario: Complex modes match the analytical closed form
- GIVEN a Rayleigh-damped SDOF or small cantilever deck with known `ζ`
- WHEN CalculixPP runs the `*FREQUENCY` + proportional `*COMPLEX FREQUENCY` job
- THEN each computed damped frequency and damping ratio SHALL match the analytical value `λ_k = -ζ_k·ω_k ± i·ω_k·√(1-ζ_k²)`, `ζ_k = (α/ω_k + β·ω_k)/2`, within the documented tolerance, and no CalculiX `*COMPLEX FREQUENCY` reference-deck comparison SHALL be asserted for this procedure

## ADDED Requirements

### Requirement: Complex-frequency keyword restrictions
The `*COMPLEX FREQUENCY` parser SHALL reject the `CORIOLIS` and `FLUTTER` keywords with an explicit "not yet implemented" error rather than silently accepting a deck it would mis-solve, because the shipped path implements only proportional-damping complex modes and does not implement the CalculiX gyroscopic Coriolis operator (which requires a rotation body load: rotor speed and axis) or the `FLUTTER` complex applied-force path.

#### Scenario: CORIOLIS keyword rejected
- GIVEN a `*COMPLEX FREQUENCY, CORIOLIS` card
- WHEN the step is parsed
- THEN the parser SHALL reject it with a clear "not yet implemented" message naming the Coriolis gyroscopic operator + rotor-speed/axis input as the missing follow-on, rather than mis-solving with proportional damping

#### Scenario: FLUTTER keyword rejected
- GIVEN a `*COMPLEX FREQUENCY, FLUTTER` card
- WHEN the step is parsed
- THEN the parser SHALL reject it with a clear "not yet implemented" message identifying the complex applied-force path as the missing follow-on
