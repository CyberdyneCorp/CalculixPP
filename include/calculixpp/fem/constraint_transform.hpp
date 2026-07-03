#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/types.hpp"

namespace cxpp::fem {

// A single free-equation contribution of a global DOF's expansion: free equation
// `eq` (0..n_free-1) weighted by `coeff`.
struct DofTerm {
  Index eq{0};
  Real coeff{0.0};
};

// The resolved expansion of one global DOF onto the master (free + prescribed) DOFs
// after eliminating multi-point-constraint slaves. A free DOF expands to itself with
// coeff 1; a slave expands to its masters' expansions with the equation coefficients;
// a prescribed (SPC) DOF contributes only to `constant` (its prescribed value, scaled
// by the chain of coefficients). The full displacement of the DOF is
//   u = constant + sum_k terms[k].coeff * u_free[terms[k].eq].
struct DofExpansion {
  std::vector<DofTerm> terms;  // free-equation contributions
  Real constant{0.0};          // prescribed contribution (SPC value carried through)
};

// The constraint transform for a model: numbers the free DOFs (SPCs and MPC slaves
// removed), and for every global DOF holds its expansion onto the free DOFs.
// Eliminating a slave DOF row/col via T^T K T (congruence) preserves symmetry and
// positive-definiteness, so the reduced SPD solve path stays usable in both the
// linear and nonlinear drivers. (spec: constraints — dependent-DOF elimination.)
//
// Over-constraint detection: constructing the transform reports a slave DOF that is
// also directly constrained (SPC) or is the dependent term of more than one equation,
// and a cyclic slave->master chain — all raise std::runtime_error with the origin(s).
struct ConstraintTransform {
  Index n_free{0};
  std::vector<Index> dof_eq;              // free equation index per DOF, or -1
  std::vector<DofExpansion> expansion;    // per global DOF (size n_dof)
  std::vector<bool> is_slave;             // per global DOF
  std::vector<std::string> messages;      // over-constraint / redundancy diagnostics

  // True if global DOF g is prescribed (SPC) — needed by callers that special-case
  // eliminated columns.
  bool prescribed(std::size_t g) const { return dof_eq[g] < 0 && !is_slave[g]; }

  // Displacement of global DOF g given the free-DOF solution `uf` and a resolved
  // per-DOF prescribed-value vector `spc` (SPC values already scaled by amplitude /
  // load factor; 0 for non-SPC DOFs). Free DOFs return uf[eq]; SPC DOFs return
  // spc[g]; MPC slaves return the equation combination of their masters. This
  // reconstructs slave motion from the current master displacements, so it is
  // correct under any load-factor scaling of the SPCs.
  Real displacement(std::size_t g, const std::vector<Real>& uf,
                    const std::vector<Real>& spc) const {
    Real v = 0.0;
    for (const DofTerm& t : expansion[g].terms) v += t.coeff * uf[static_cast<std::size_t>(t.eq)];
    // Add the prescribed (SPC master) contribution, recomputed from `spc` so it
    // tracks the current load factor (expansion.constant used the full magnitude).
    v += spc_contribution(g, spc);
    return v;
  }

 private:
  // The SPC (prescribed-master) part of DOF g's expansion, using the supplied
  // current SPC values. For a free DOF it is 0 (no SPC master); for an SPC DOF it is
  // spc[g]; for a slave it is the coefficient-weighted sum over its expansion's SPC
  // masters. Computed by re-walking is not needed: expansion.constant already stores
  // this with the FULL prescribed magnitude, so we scale per-DOF via `spc` only when
  // the caller needs load-factor tracking. Here we recompute from spc via the stored
  // per-DOF SPC weights.
 public:
  // Per global DOF: the SPC masters (global dof, weight) that contribute to its
  // displacement. Filled by build_constraint_transform alongside `expansion`.
  std::vector<std::vector<DofTerm>> spc_terms;  // uses DofTerm{eq=global spc dof, coeff}

  Real spc_contribution(std::size_t g, const std::vector<Real>& spc) const {
    Real v = 0.0;
    for (const DofTerm& t : spc_terms[g]) v += t.coeff * spc[static_cast<std::size_t>(t.eq)];
    return v;
  }
};

// Build the constraint transform for a model. `spc_constrained` marks the DOFs pinned
// by *BOUNDARY (with values in `prescribed_value`). The equations come from
// model.expand_constraints(). Throws std::runtime_error on an over-constraint
// (slave also an SPC, doubly-defined dependent DOF, or a dependency cycle).
ConstraintTransform build_constraint_transform(
    const Model& model, const std::vector<bool>& spc_constrained,
    const std::vector<Real>& prescribed_value);

}  // namespace cxpp::fem
