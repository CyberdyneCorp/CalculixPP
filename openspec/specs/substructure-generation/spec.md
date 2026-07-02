# Substructure Generation and Matrix Export

## Purpose

CalculixPP generates reduced substructures (superelements) and exports system
matrices, ported to pure C++20. A `*SUBSTRUCTURE GENERATE` step condenses the model
onto a set of retained (master) DOFs via static (Guyan) or Craig-Bampton reduction,
producing reduced stiffness and mass operators for reuse. The condensation and any
partial factorization run through **NumPP** (dense/sparse reduction, not SPOOLES/
PARDISO), and assembly/solve hot paths route through the `ComputeBackend` (CPU/NumPP by
default, GPU optional and never required). Exported reduced and global matrices SHALL
match reference CalculiX within tolerance and are reachable from the Python bindings.
(ref: src/CalculiX.c substructure path, src/substructures.f)

**Porting Phase:** 4 — Dynamics & eigenproblems

## Requirements

### Requirement: Substructure generation over retained DOFs
A `*SUBSTRUCTURE GENERATE` step SHALL condense the model onto the retained DOFs declared by `*RETAINED NODAL DOFS`, producing a superelement whose reduced operators expose only those boundary DOFs. The step SHALL support static (Guyan) reduction and, when a mass matrix is present, Craig-Bampton reduction retaining fixed-interface modes, and SHALL be reachable from the Python bindings. (ref: src/substructures.f)

#### Scenario: Generate a superelement
- GIVEN a model with a `*SUBSTRUCTURE GENERATE` step and `*RETAINED NODAL DOFS`
- WHEN the step runs
- THEN the solver SHALL emit a superelement whose degrees of freedom are exactly the retained DOFs

#### Scenario: Craig-Bampton with mass
- GIVEN a `*SUBSTRUCTURE GENERATE` step requesting a mass matrix and fixed-interface modes
- WHEN the step runs
- THEN the reduced model SHALL include the retained-DOF condensation plus the requested fixed-interface modal DOFs

### Requirement: Retained nodal DOFs definition
The `*RETAINED NODAL DOFS` card SHALL define the set of boundary (master) DOFs kept in the reduction, given as node or node-set entries with first and last DOF ranges, and every retained DOF SHALL survive as an active DOF of the exported superelement. (ref: src/substructures.f)

#### Scenario: DOF range per node set
- GIVEN a `*RETAINED NODAL DOFS` card listing a node set with DOFs 1 through 3
- WHEN the substructure is generated
- THEN translational DOFs 1, 2, and 3 of each node in the set SHALL be retained

#### Scenario: No retained DOFs declared
- GIVEN a `*SUBSTRUCTURE GENERATE` step with no `*RETAINED NODAL DOFS` in the step
- WHEN the job runs
- THEN the program SHALL report an error that the retained DOF set is empty

### Requirement: Reduced matrix reduction via NumPP
Static condensation of the interior DOFs SHALL be performed with NumPP, factoring the interior-interior block and forming the Schur-complement reduced stiffness (and mass, when requested) on the retained DOFs, executed through the `ComputeBackend`. (ref: src/substructures.f)

#### Scenario: No GPU present
- GIVEN a build with no CUDA/OpenCL/Metal backend available
- WHEN a `*SUBSTRUCTURE GENERATE` step runs
- THEN the interior condensation SHALL be factored and reduced on the CPU/NumPP backend and produce correct reduced operators

### Requirement: Reduced and global matrix export
The `*SUBSTRUCTURE MATRIX OUTPUT` card SHALL write the reduced stiffness and mass matrices of the superelement for reuse, and `*MATRIX ASSEMBLE` SHALL assemble and export the global system stiffness and mass matrices to file. (ref: src/substructures.f, src/CalculiX.c substructure path)

#### Scenario: Write reduced matrices
- GIVEN a `*SUBSTRUCTURE MATRIX OUTPUT` card requesting stiffness and mass
- WHEN the substructure step completes
- THEN the reduced stiffness and mass matrices SHALL be written to the substructure output file ordered by retained DOF

#### Scenario: Export global matrices
- GIVEN a step with `*MATRIX ASSEMBLE` requesting stiffness and mass
- WHEN the step runs
- THEN the assembled global stiffness and mass matrices SHALL be exported to file (see linear-algebra-and-solvers for the assembled operator representation)

### Requirement: Reference-result fidelity
For every substructure `test/` deck, the exported reduced and global matrices SHALL match reference CalculiX output within a documented relative tolerance, validated through the pytest bindings.

#### Scenario: Reduced matrix matches reference
- GIVEN a reference CalculiX substructure deck and its expected reduced matrices
- WHEN CalculixPP runs the same deck
- THEN each exported reduced matrix entry SHALL match the reference within the documented relative tolerance
