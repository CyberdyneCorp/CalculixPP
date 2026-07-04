# Input Deck Parsing

## Purpose

CalculixPP reads the same Abaqus-style ASCII input deck as reference CalculiX,
tokenizes its keyword cards, and builds the in-memory model and history the solver
operates on. This capability defines the lexical/structural rules of the deck and the
dispatch of keyword cards; parsed results SHALL match reference CalculiX for the
corresponding `test/` deck. Implemented as a portable C++20 tokenizer/parser with no
Fortran fixed-field assumptions. (ref: src/CalculiX.c, src/calinput.f, src/keystart.f)

**Porting Phase:** 1 — Foundation & vertical slice (parses the linear-static keyword subset; later phases extend the recognized card set).
## Requirements
### Requirement: Job invocation from an input file
CalculixPP SHALL accept a job name (`-i jobname` on the CLI, or an equivalent argument
from the Python bindings) and read the model from `jobname.inp`, producing result
files prefixed with the same job name. The same entry point SHALL be reachable from the
Python bindings (see python-bindings spec).

#### Scenario: Run a job
- GIVEN an input file `model.inp`
- WHEN the user runs the solver with job name `model`
- THEN the program SHALL parse `model.inp`
- AND write results to `model.frd` and `model.dat`

### Requirement: Keyword-card structure
The deck SHALL be organized as keyword cards, each beginning with a line whose first
non-blank character is `*` (a card name), optionally followed by comma-separated
parameters, and followed by data lines until the next card. Lines beginning with `**`
SHALL be treated as comments and ignored. (ref: src/calinput.f, src/keystart.f)

#### Scenario: Comment line
- GIVEN a line starting with `**`
- WHEN the parser reads it
- THEN the line SHALL be ignored and not interpreted as data

#### Scenario: Unknown keyword
- GIVEN a card name that is not recognized
- WHEN the parser encounters it
- THEN the program SHALL report an error identifying the unknown keyword and its line

### Requirement: Model data vs. history data separation
The deck SHALL consist of model data (geometry, materials, sections, sets) followed by
one or more analysis steps. Each analysis step SHALL be delimited by a `*STEP` card and
a `*END STEP` card; history cards (procedure, loads, output requests) SHALL appear
inside a step. (ref: src/calinput.f, src/steps.f)

#### Scenario: History card outside a step
- GIVEN a procedure card such as `*STATIC` placed outside any `*STEP`/`*END STEP` pair
- WHEN the parser processes it
- THEN the program SHALL report an error

### Requirement: File inclusion
The parser SHALL support `*INCLUDE, INPUT=file` to splice the contents of another file
into the deck at the point of the card, resolving nested includes. (ref: src/calinput.f)

#### Scenario: Included file
- GIVEN `*INCLUDE, INPUT=mesh.inp` in the main deck
- WHEN the parser reaches the card
- THEN the contents of `mesh.inp` SHALL be parsed as if written inline at that point

### Requirement: Case-insensitivity and free-field input
Keyword names and parameters SHALL be interpreted case-insensitively, and numeric data
SHALL be read in free (comma-separated) field format independent of column position.

#### Scenario: Mixed case keyword
- GIVEN `*Static` and `*STATIC`
- WHEN parsed
- THEN both SHALL select the static procedure

### Requirement: Deterministic parse against reference
For a given `test/` deck, the parsed model and history SHALL be equivalent to the
structures reference CalculiX builds from the same deck, so downstream results match
within the documented numerical tolerance.

#### Scenario: Reference-equivalent parse
- GIVEN a `test/` deck parsed by both reference CalculiX and CalculixPP
- WHEN the resulting node, element, set, and step data are compared
- THEN they SHALL be equivalent up to ordering and floating-point tolerance

### Requirement: Model-control cards (Phase 1)
The parser SHALL support `*HEADING` to capture a free-text model title/heading for the job, and `*NO ANALYSIS` (card `*NOANALYSIS`) to parse and check the deck without performing a solve, matching reference CalculiX semantics. (ref: src/headings.f, src/noanalysiss.f, src/calinput.f)

#### Scenario: Free-text heading
- GIVEN a `*HEADING` card followed by one or more text lines
- WHEN the parser reads the card
- THEN the text SHALL be stored as the model heading and carried to the output files

#### Scenario: Parse-only run
- GIVEN a deck containing `*NO ANALYSIS`
- WHEN the job runs
- THEN the deck SHALL be fully parsed and checked
- AND no solve SHALL be performed

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

