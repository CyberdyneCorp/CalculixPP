#pragma once
#include <array>
#include <vector>

#include "calculixpp/core/types.hpp"

namespace cxpp {

using Voigt6 = std::array<Real, 6>;  // stress/strain: xx, yy, zz, xy, xz, yz

// Nodal result fields for a static step, aligned with mesh node indices.
struct StaticFields {
  std::vector<Vec3> displacement;  // U
  std::vector<Voigt6> stress;      // S (averaged nodal stress)
  std::vector<Voigt6> strain;      // E (averaged nodal strain, engineering shear)
  std::vector<Vec3> reaction;      // RF (f_int - f_ext)
};

// Nodal result fields for a heat-transfer step, aligned with mesh node indices
// (spec: heat-transfer-analysis / results-output). `temperature` is the NT field;
// `flux_reaction` is the concentrated nodal heat-flux reaction RFL at prescribed-
// temperature nodes (Kt*T at the constrained DOFs), analogous to mechanical RF.
struct ThermalFields {
  std::vector<Real> temperature;     // NT
  std::vector<Real> flux_reaction;   // RFL
};

// Combined result of a *COUPLED TEMPERATURE-DISPLACEMENT step (spec:
// heat-transfer-analysis — coupled). `thermal` is the solved temperature field
// (NT/RFL); `mechanical` is the displacement/stress/reaction field computed with the
// thermal strain eps_th = alpha (T - Tref) applied. The one-way (sequential) scheme
// fills `thermal` first, then `mechanical` from the resulting temperature field.
struct CoupledFields {
  ThermalFields thermal;
  StaticFields mechanical;
};

}  // namespace cxpp
