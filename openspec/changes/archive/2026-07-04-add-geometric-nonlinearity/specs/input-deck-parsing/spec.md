## ADDED Requirements

### Requirement: *BUCKLE procedure card
The parser SHALL accept a `*BUCKLE` card, set the buckling procedure, and read the requested number of buckling modes from its data line (accepting and ignoring the ARPACK-style tolerance/subspace/iteration fields), instead of rejecting it as a deferred capability.

The card sets `Procedure::Buckling` and records `num_buckling_modes` from the first data field, as the
`*FREQUENCY` card records `num_eigenvalues`. The trailing accuracy/subspace/iteration fields are
recognized but do not affect the model. (ref: src/bucklings.f, src/calinput.f)

#### Scenario: Records the buckling procedure and mode count
- GIVEN a deck containing `*BUCKLE` with a mode count on its data line
- WHEN parsed
- THEN the model procedure SHALL be `Buckling`
- AND the requested mode count SHALL be recorded

#### Scenario: Trailing ARPACK-style fields are ignored
- GIVEN `*BUCKLE` with trailing tolerance/subspace/iteration fields (e.g. `10, 1.e-2`)
- WHEN parsed
- THEN those fields SHALL be accepted and ignored without error
