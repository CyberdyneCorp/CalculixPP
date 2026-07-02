# Results and Output

## Purpose

Writes CalculixPP computed results to files for postprocessing and inspection: nodal
and element field results to the `.frd` file (CalculiX/CGX-compatible format), tabular
results to the `.dat` file, and supports restart. Output requests select which
quantities are written and where. In Phase 1 the `.frd`/`.dat` writers are the
end-of-pipeline stage that persists the linear-static solution; the same results are
also reachable through the Python bindings. Format and label semantics mirror the
reference; only the implementation stack changes. (ref: src/frd.c, ref: src/frdheader.c,
ref: src/printout.f, ref: src/restartwrite.f)

**Porting Phase:** `.frd`/`.dat` writing for the linear-static vertical slice
(displacements U, stresses S, strains E, reaction forces RF) = Phase 1; binary `.frd`,
restart, section prints, and thermal/field labels beyond the Phase-1 set = later phases.

## Requirements

### Requirement: FRD result file
CalculixPP SHALL write nodal and element results to `jobname.frd` in the CalculiX
results format readable by the CGX postprocessor. The quantities written SHALL be
controlled by `*NODE FILE`, `*EL FILE`, `*NODE OUTPUT`, and `*ELEMENT OUTPUT` requests
with selectable result labels; the Phase-1 label set SHALL include at least `U`
(displacements), `S` (stresses), `E` (strains), and `RF` (reaction forces). (ref: src/frd.c)

#### Scenario: Request displacements and stresses
- GIVEN `*NODE FILE` with `U` and `*EL FILE` with `S` in a step
- WHEN the step completes an output increment
- THEN displacements (U) and stresses (S) SHALL be written to `jobname.frd` for that increment
- AND the file SHALL be readable by the CGX postprocessor

### Requirement: Linear-static results output
For a completed linear-static solve, CalculixPP SHALL write the recovered nodal
displacements and element stresses/strains to `jobname.frd` and any requested tabular
quantities to `jobname.dat`, closing the Phase-1 pipeline. The written results SHALL be
consistent with the same results exposed through the Python bindings. (ref: src/frd.c,
ref: src/printout.f)

#### Scenario: Persist the vertical-slice solution
- GIVEN a linear `*STATIC` step that has solved `K u = f` and recovered element results
- WHEN the output increment is written
- THEN nodal displacements SHALL be written to `jobname.frd`
- AND element stresses and strains SHALL be written to `jobname.frd`
- AND requested tabular values SHALL be written to `jobname.dat`
- AND the same displacement and stress data SHALL be retrievable via the Python bindings

### Requirement: DAT print file
CalculixPP SHALL write requested tabular results to `jobname.dat` via `*NODE PRINT`,
`*EL PRINT`, and `*SECTION PRINT`, including totals and section results where requested.
(ref: src/printout.f, ref: src/printoutface.f)

#### Scenario: Nodal print
- GIVEN a `*NODE PRINT` request for `U`
- WHEN the step completes an output increment
- THEN the nodal displacement table SHALL be written to `jobname.dat`

#### Scenario: Section print
- GIVEN a `*SECTION PRINT` request on a surface for `SOF` (section forces)
- WHEN the step runs
- THEN the integrated section forces SHALL be written to `jobname.dat`

### Requirement: Output frequency control
Output requests SHALL honor a frequency parameter (`FREQUENCY=`) and `*TIME POINTS`
so results are written at selected increments or times rather than every increment.
(ref: src/calinput.f, ref: src/timepointss.f)

#### Scenario: Write every Nth increment
- GIVEN `*NODE FILE, FREQUENCY=5`
- WHEN the step runs
- THEN results SHALL be written only on every 5th increment

### Requirement: Restart
CalculixPP SHALL write restart data via `*RESTART, WRITE` and resume a job from saved
state via `*RESTART, READ`, preserving the model and solution state across runs.
(ref: src/restartwrite.f, ref: src/restartread.f, ref: src/restarts.f)

#### Scenario: Resume from restart
- GIVEN a previous run wrote restart data at the end of a step
- WHEN a new job specifies `*RESTART, READ`
- THEN the new job SHALL continue from the saved state rather than starting from scratch

### Requirement: Output formats
The `.frd` output SHALL support both ASCII and binary formats, selectable via the
`OUTPUT` parameter on `*NODE FILE`/`*EL FILE` (or the global output setting). ASCII
output SHALL be available in Phase 1 as the CGX-compatible default. (ref: src/frd.c,
ref: src/frdheader.c)

#### Scenario: Binary FRD
- GIVEN an output request selecting binary format
- WHEN results are written
- THEN the `.frd` file SHALL be written in binary rather than ASCII

### Requirement: Global output control (Phase 1)
CalculixPP SHALL support `*OUTPUT` to set global output options, selecting `.frd` versus `.dat` output and controlling the 2-D-to-3-D expansion of results for plane/axisymmetric/shell/beam elements, so requested quantities are written in the expanded 3-D representation where applicable, matching reference CalculiX semantics. (ref: src/frd.c, ref: src/calinput.f)

#### Scenario: Expand 2-D results to 3-D
- GIVEN a plane-stress model with `*OUTPUT` requesting 3-D expansion
- WHEN results are written
- THEN the `.frd` results SHALL be written in the expanded 3-D representation of the elements

#### Scenario: Select global output target
- GIVEN a `*OUTPUT` card selecting the output file target
- WHEN a subsequent output request is written
- THEN the results SHALL be directed to the selected target per the global setting
