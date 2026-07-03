#pragma once
#include <string>
#include <vector>

#include "calculixpp/core/types.hpp"

namespace cxpp {

// One term of a linear multi-point constraint: coefficient * u_(node_id, comp).
// `comp` is 1..3 (structural DOF). (spec: constraints — *EQUATION.)
struct EquationTerm {
  Index node_id{};
  int comp{};      // 1..3
  Real coeff{0.0};
};

// A linear homogeneous constraint  sum_i coeff_i * u_(node_i, comp_i) = 0.
// The FIRST term is the dependent (slave) DOF, eliminated during assembly via an
// SPD-preserving congruence transform (u_slave expressed from the master terms).
// (spec: constraints — Linear equations / *EQUATION, *MPC, *RIGID BODY, *COUPLING,
// *TIE all reduce to a list of these.)
struct Equation {
  std::vector<EquationTerm> terms;  // terms[0] is the dependent DOF
  std::string origin;               // human-readable source (e.g. "*EQUATION", "*MPC BEAM")
};

}  // namespace cxpp
