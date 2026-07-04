#pragma once
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/types.hpp"
#include "calculixpp/fem/constraint_transform.hpp"
#include "calculixpp/fem/material_model.hpp"

namespace cxpp::fem {

// Assembled reduced linear system for a linear-static step.
// Holds the free/free stiffness in COO form (duplicates already summed) and the
// reduced load vector, plus the DOF map needed to expand the solution back to
// nodal displacements. (spec: static-analysis — linear static end-to-end pipeline.)
struct LinearSystem {
  Index n_free{0};
  std::vector<Index> rows;  // COO row indices (0..n_free-1)
  std::vector<Index> cols;  // COO col indices
  std::vector<Real> vals;   // COO values
  std::vector<Real> rhs;    // reduced load vector, size n_free

  // dof_eq[node_index*3 + comp] = free equation (>=0) or -1 if constrained (SPC or
  // MPC slave). Mirrors transform.dof_eq for callers that only need the free numbering.
  std::vector<Index> dof_eq;
  // prescribed[node_index*3 + comp] = prescribed value for SPC DOFs (else 0).
  std::vector<Real> prescribed;
  // prescribed_amp[node_index*3 + comp] = *AMPLITUDE name for that constrained DOF
  // (empty -> default linear ramp). Lets the driver scale prescribed BCs over time.
  std::vector<std::string> prescribed_amp;

  // Multi-point-constraint transform: expands each global DOF onto the free DOFs
  // (SPCs + MPC slaves eliminated). Empty terms / identity when the model has no
  // constraints. Diagnostics (over-constraint / redundancy) live in transform.messages.
  ConstraintTransform transform;
};

// Assemble K_ff u_f = f_ext,f - K_fc u_c and the reduced load, using
// static condensation of the single-point constraints.
LinearSystem assemble_linear_static(const Model& model);

// Assemble the global mass matrix on the SAME free-DOF numbering as
// assemble_linear_static (spec: modal-and-buckling-analysis — *FREQUENCY). Each
// element's consistent (or lumped) mass M_e = ∫rho Nᵀ N dV is scattered through the
// constraint transform exactly like the stiffness, so the returned COO triplets and
// the free-DOF count match the stiffness system produced by assemble_linear_static
// for the SAME model — the generalized eigenproblem K x = λ M x is well-posed on the
// shared numbering. `*MASS` point-mass connectors add a diagonal nodal mass on their
// three translational DOFs. The returned LinearSystem carries the same dof_eq /
// prescribed / transform data as the stiffness system; rhs is unused (left zero).
// `lumped` selects row-sum lumped element mass instead of the consistent matrix.
// (Deactivated elements — *MODEL CHANGE, REMOVE — carry no mass, mirroring K.)
LinearSystem assemble_mass(const Model& model, bool lumped = false);

// Per-element material models and their per-integration-point state, aligned with
// mesh.elements(). Built once by make_material_points and advanced across Newton
// iterations by the nonlinear driver. `state[e]` has one MaterialState per Gauss
// point of element e. (spec: material-models — material-point state carried across
// increments.)
struct MaterialPoints {
  std::vector<std::unique_ptr<MaterialModel>> models;
  std::vector<std::vector<MaterialState>> state;

  // Commit every integration point's trial history into its committed slot. Called
  // by the nonlinear driver once an increment converges, so the next increment /
  // Newton iteration return-maps from the converged state (spec: material-models —
  // trial state kept per iteration, committed on increment convergence).
  void commit() {
    for (auto& elem_state : state)
      for (MaterialState& s : elem_state) s.commit();
  }
};

// Build the material-point store for a model: one MaterialModel + a Gauss-point
// state vector per element. Currently every element is linear-elastic
// (ElasticIsoMaterial from element_elastic()).
MaterialPoints make_material_points(const Model& model);

// Material-point tangent assembly at the current full nodal displacement `u`
// (size mesh.num_nodes()). For each element it forms B ue, evaluates the material at
// each Gauss point (advancing mp.state), and accumulates the consistent tangent into
// the reduced free/free COO system and the internal force into `f_int` (size
// 3*num_nodes). The returned LinearSystem carries the same DOF map / prescribed data
// as assemble_linear_static and an rhs holding the reduced external load at full
// magnitude minus the eliminated tangent columns; the Newton driver overwrites rhs
// with its residual before solving. For linear elasticity at u=0 the tangent equals
// assemble_linear_static's K exactly. (spec: nonlinear-solution-control — tangent +
// residual from material-point assembly.)
LinearSystem assemble_material_tangent(const Model& model,
                                       const std::vector<Vec3>& u,
                                       MaterialPoints& mp,
                                       std::vector<Real>& f_int);

}  // namespace cxpp::fem
