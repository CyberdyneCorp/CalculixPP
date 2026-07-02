#pragma once
#include "calculixpp/core/model.hpp"
#include "calculixpp/core/results.hpp"

namespace cxpp::fem {

// Given a solved displacement field (fields.displacement, all nodes), recover the
// averaged nodal stresses and the reaction forces. Integration-point stresses are
// extrapolated to element nodes (C3D4: constant; C3D10: a4 corner extrapolation +
// midside averaging) and averaged across elements. RF = f_int - f_ext.
// (spec: static-analysis / results-output — reimplemented from CalculiX
// resultsmech.f / extrapolate.f, not copied.)
void recover_fields(const Model& model, StaticFields& fields);

}  // namespace cxpp::fem
