# Adaptive Mesh Refinement

## Purpose

Refines an existing tetrahedral mesh during an analysis so regions with high solution
error (or a chosen field variable) are meshed more finely, then re-runs the analysis on
the refined mesh. This is solve-driven refinement of a mesh that already exists in the
model — distinct from the `mesh-processing` capability, which generates a mesh from CAD
geometry. In CalculixPP the tetrahedral create/split/sliver-removal operations are
performed through **CyberCadKernel** remeshing primitives; the driving error field comes
from a completed solve. Results on the refined mesh SHALL match reference CalculiX within
a documented tolerance. (ref: src/refinemesh.c, src/refinemeshs.f, src/generatetet_refine.f)

**Porting Phase:** 5 — Advanced physics & adaptivity (cross-cuts static/heat-transfer once those phases land).

## Requirements

### Requirement: Refine-mesh request
CalculixPP SHALL support `*REFINE MESH` inside a step to request refinement of the current tetrahedral mesh (linear C3D4 or quadratic C3D10), driven by a selected field label such as an error indicator or result field. (ref: src/refinemeshs.f)

#### Scenario: Request refinement on a field
- GIVEN a step containing `*REFINE MESH` referencing a field label
- WHEN the step completes
- THEN CalculixPP SHALL build a new mesh refined where that field is large and store it for a subsequent solve

#### Scenario: Non-tetrahedral mesh
- GIVEN a model whose elements are not tetrahedra
- WHEN `*REFINE MESH` is requested
- THEN CalculixPP SHALL report that mesh refinement requires a tetrahedral mesh

### Requirement: Solve–refine–resolve loop
When refinement is requested, CalculixPP SHALL, after creating the refined mesh, restart the calculation from the beginning on the new mesh so the analysis is ultimately solved on the refined mesh. (ref: src/CalculiX.c refinement loop)

#### Scenario: Automatic restart on the refined mesh
- GIVEN a completed step that produced a refined mesh
- WHEN the refinement loop advances
- THEN the model data SHALL be reloaded and the analysis SHALL re-run on the refined mesh

### Requirement: Tetrahedral mesh quality
The refinement process SHALL maintain a valid tetrahedral mesh — creating and removing tetrahedra and removing degenerate (sliver) elements via CyberCadKernel — so the resulting mesh is usable by the solver, reusing the quality metrics defined in the mesh-processing capability. (ref: src/createtet.f, src/cattet.f, src/removetet_sliver.f)

#### Scenario: Sliver removal
- GIVEN refinement produces poorly shaped (sliver) tetrahedra
- WHEN the mesh is finalized
- THEN such tetrahedra SHALL be removed or repaired before the mesh is written

### Requirement: Base and refined mesh selection
CalculixPP SHALL support `*INITIAL MESH` to mark the base mesh and `*USE REFINED MESH` to run the analysis on a previously produced refined mesh. (ref: src/calinput.f)

#### Scenario: Reuse a previous refined mesh
- GIVEN a `*USE REFINED MESH` card and a previously generated refined mesh
- WHEN the job runs
- THEN the analysis SHALL use the refined mesh rather than the original base mesh

### Requirement: Refined mesh output
CalculixPP SHALL write the newly generated mesh (nodes and elements, including quadratic mid-side nodes when the source mesh is quadratic) so it can be reused or inspected, and this write SHALL be reachable from the Python bindings. (ref: src/writenewmesh.c, src/writerefinemesh.f)

#### Scenario: Preserve element order
- GIVEN refinement of a quadratic (C3D10) mesh
- WHEN the new mesh is written
- THEN the written elements SHALL include mid-side nodes so the refined mesh stays quadratic

### Requirement: Reference-result fidelity
The solution obtained on a refined mesh for a given reference `test/` deck SHALL match the reference CalculiX result within the documented numerical tolerance, on the CPU backend with no GPU present. (ref: src/refinemesh.c)

#### Scenario: Refined-mesh result agreement
- GIVEN a reference `*REFINE MESH` deck and its reference CalculiX result
- WHEN CalculixPP runs the solve–refine–resolve loop on the CPU backend
- THEN the final displacements/stresses SHALL agree with the reference within the documented tolerance
