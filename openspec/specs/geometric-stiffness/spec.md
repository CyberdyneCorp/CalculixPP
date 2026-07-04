# geometric-stiffness Specification

## Purpose
TBD - created by archiving change add-geometric-nonlinearity. Update Purpose after archive.
## Requirements
### Requirement: Initial-stress (geometric) stiffness operator
CalculixPP SHALL assemble a geometric (initial-stress) stiffness matrix `K_geo` about a reference configuration from a supplied per-integration-point reference stress field, on the same free-DOF numbering and constraint transform as the linear stiffness, where each element contributes the block-diagonal initial-stress term `K_geo[3a+i][3b+j] = δ_ij · Σ_gp (∇N_a · σ · ∇N_b) · det(J) · w`.

The operator uses the SAME Gauss rule as the element stiffness for consistency, is symmetric, and
introduces no hourglass term (only the stress-gradient product). This is the `buckling=1` branch of
reference CalculiX `e_c3d.f`, which places `K_geo` into the pencil with a negative sign. (ref: src/e_c3d.f, src/bucklings.f)

#### Scenario: Shares the free-DOF layout of the linear stiffness
- GIVEN a model and a per-element per-Gauss reference stress field
- WHEN `K_geo` is assembled
- THEN it SHALL be symmetric
- AND it SHALL share the free-DOF count and `dof_eq` of `assemble_linear_static` for the same model

#### Scenario: Zero stress field yields a zero matrix
- GIVEN a zero reference stress field
- WHEN `K_geo` is assembled
- THEN every entry SHALL be zero
- AND the linear-static path SHALL be unaffected

#### Scenario: Single element matches the analytic initial-stress matrix
- GIVEN a single element under a known uniform stress
- WHEN its element `K_geo` is formed
- THEN it SHALL equal the analytic block-diagonal `∇N·σ·∇N` initial-stress matrix

### Requirement: Two-step prestress driver
CalculixPP SHALL compute a geometric-stiffness eigenbasis by first solving the linear static response to a reference load, recovering the reference integration-point stresses, then assembling `K` and `K_geo` and solving the generalized eigenproblem, without altering the linear-static path for non-buckling decks.

The Cholesky anchor for the generalized reduction is the UNLOADED `K` (SPD before the reference load
is applied), the classical linear-buckling formulation that keeps the reduction well-posed while
`−K_geo` is indefinite. Reference stress recovery reuses the shipped stress path without changing its
existing output. (ref: src/arpackbu.c, src/resultsini.f)

#### Scenario: Prestress solve feeds the geometric stiffness
- GIVEN a `*BUCKLE` step with a reference load
- WHEN it runs
- THEN a linear prestress solve SHALL produce the reference stress field feeding `K_geo`

#### Scenario: Non-buckling decks are untouched
- GIVEN a deck with no `*BUCKLE` / geometric step
- WHEN it is solved
- THEN the linear-static and Newton paths SHALL be byte-for-byte unchanged
- AND the existing stress-recovery output SHALL be byte-identical

