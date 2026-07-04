#pragma once
#include <cstddef>
#include <vector>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/types.hpp"
#include "calculixpp/numerics/eigensolution.hpp"

// Linear-buckling two-step prestress driver (spec: geometric-stiffness — two-step
// prestress driver; modal-and-buckling-analysis — *BUCKLE). Solves the buckling pencil
// (K + λ K_geo) φ = 0 for the lowest positive load factors:
//   Step A — solve the linear static response K u0 = f_ref to the step's reference
//            load, recover the reference Gauss-point stress field σ_ref(u0).
//   Step B — assemble K and K_geo(σ_ref) and extract the lowest buckling factors
//            (ascending positive) + their mode shapes via the dense generalized path.
// The static/Newton drivers are untouched: this reuses solve_linear_static +
// recover_gauss_stress + assemble_geometric_stiffness + extract_buckling_modes.
namespace cxpp::numerics {

// Result of a *BUCKLE step: the lowest buckling load factors (ascending positive) and
// the full eigenbasis (mode shapes, free-DOF numbering) for output plumbing. The
// critical load is factors[0] * f_ref.
struct BucklingReport {
  std::vector<Real> factors;  // buckling load factors λ_i (ascending positive)
  EigenBasis basis;           // mode shapes + free-DOF numbering (eigenvalue = λ)
};

// Run the two-step buckling analysis for `model`, extracting the lowest `num_modes`
// positive factors. Uses the dense generalized eigen path (validatable at beamb scale).
BucklingReport solve_buckling(const Model& model, std::size_t num_modes);

}  // namespace cxpp::numerics
