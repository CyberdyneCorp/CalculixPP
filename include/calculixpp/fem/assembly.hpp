#pragma once
#include <string>
#include <vector>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/types.hpp"

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

  // dof_eq[node_index*3 + comp] = free equation (>=0) or -1 if constrained.
  std::vector<Index> dof_eq;
  // prescribed[node_index*3 + comp] = prescribed value for constrained DOFs (else 0).
  std::vector<Real> prescribed;
  // prescribed_amp[node_index*3 + comp] = *AMPLITUDE name for that constrained DOF
  // (empty -> default linear ramp). Lets the driver scale prescribed BCs over time.
  std::vector<std::string> prescribed_amp;
};

// Assemble K_ff u_f = f_ext,f - K_fc u_c and the reduced load, using
// static condensation of the single-point constraints.
LinearSystem assemble_linear_static(const Model& model);

}  // namespace cxpp::fem
