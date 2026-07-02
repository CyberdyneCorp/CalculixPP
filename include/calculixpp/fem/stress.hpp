#pragma once
#include <vector>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/results.hpp"
#include "calculixpp/core/types.hpp"

namespace cxpp::fem {

// Global internal force vector f_int = sum_e Ke_e u_e (size 3*num_nodes, DOF order
// [u,v,w] per node). For linear elasticity Ke_e is the element stiffness, so this is
// the exact internal force that the Newton residual r = lambda*f_ext - f_int uses;
// it is structured to accept a nonlinear tangent/stress later. `u` is the current
// nodal displacement (all nodes). (spec: nonlinear-solution-control — residual.)
std::vector<Real> internal_force(const Model& model, const std::vector<Vec3>& u);

// Given a solved displacement field (fields.displacement, all nodes), recover the
// averaged nodal stresses and the reaction forces. Integration-point stresses are
// extrapolated to element nodes (C3D4: constant; C3D10: a4 corner extrapolation +
// midside averaging) and averaged across elements. RF = f_int - f_ext.
// (spec: static-analysis / results-output — reimplemented from CalculiX
// resultsmech.f / extrapolate.f, not copied.)
void recover_fields(const Model& model, StaticFields& fields);

}  // namespace cxpp::fem
