## MODIFIED Requirements

### Requirement: Linear buckling analysis
A `*BUCKLE` step SHALL compute buckling load factors and buckling mode shapes by solving `(K + λ K_geo) φ = 0`, where `K_geo` derives from a reference load's prestress stress state, returning the buckling factors in ascending positive order with their mode shapes, reachable from the CLI and Python bindings.

The eigenproblem is solved as the pencil `(A = −K_geo, B = K)` with `B = K` SPD (the unloaded
stiffness), mapping `λ = −1/θ` and filtering to positive factors; the critical load is the smallest
positive factor times the reference load. The dense generalized path (NumPP Cholesky of `K` + `eigh`)
is the primary validatable path; the scalable sparse path is gated on upstream SciPP target-selection
support. (ref: src/arpackbu.c, src/bucklings.f, src/e_c3d.f)

#### Scenario: Buckling factors and mode shapes
- GIVEN a `*BUCKLE` step with a reference load
- WHEN the step runs
- THEN the solver SHALL output buckling load factors in ascending positive order with their mode shapes
- AND the critical load SHALL be the smallest positive factor times the reference load

#### Scenario: Matches the beamb reference deck
- GIVEN the stock `beamb` buckling deck (C3D20R Euler column, `*BUCKLE 10, 1.e-2`)
- WHEN CalculixPP runs it
- THEN the reported buckling factors SHALL match `beamb.dat.ref` (λ_1 = 48.15, λ_2 = 106.3) within a documented relative tolerance

#### Scenario: Matches the Euler critical load
- GIVEN a slender column meshed with solid elements under axial compression
- WHEN a `*BUCKLE` step runs
- THEN the lowest factor times the reference load SHALL match the Euler critical load `π²EI/(kL)²` within tolerance for pinned-pinned (`k=1`) and clamped-free (`k=2`)
