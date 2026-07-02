#pragma once
#include <array>
#include <vector>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/types.hpp"

namespace cxpp::fem {

// Local element node indices (0-based) of a tetrahedral face (1..4), CalculiX
// ifacet ordering: 3 corners for C3D4, 3 corners + 3 midsides for C3D10.
std::vector<int> face_nodes(ElementType type, int face);

// Global external load vector (size 3*num_nodes) at full magnitude: concentrated
// loads (*CLOAD), consistent nodal forces from pressure (*DLOAD P<face> / *DSLOAD),
// and body loads (*DLOAD GRAV / CENTRIF). Amplitude references are evaluated at the
// end of the step (factor 1), so this is the fully-applied load.
std::vector<Real> external_load_vector(const Model& model);

// Same, but scaled per load by the step fraction `lambda` in [0,1]: each load is
// multiplied by amplitude_factor(load.amplitude, lambda * total), which reduces to
// `lambda` for loads without an amplitude. Used by the nonlinear driver so
// time-varying (*AMPLITUDE) loads are applied at the correct step fraction, while
// unamplified loads keep the linear 0->1 ramp. (spec: nonlinear-solution-control —
// Amplitude-driven time stepping.)
std::vector<Real> external_load_vector(const Model& model, Real lambda);

}  // namespace cxpp::fem
