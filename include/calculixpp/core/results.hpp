#pragma once
#include <array>
#include <vector>

#include "calculixpp/core/types.hpp"

namespace cxpp {

using Voigt6 = std::array<Real, 6>;  // stress/strain: xx, yy, zz, xy, xz, yz

// Contact result at one slave node (spec: contact — contact output CSTR). CalculiX writes
// contact stress/status per slave node; this carries the node id, the OPEN/CLOSED status
// (closed when the node is in contact / penetrating), the normal contact pressure `p`
// (>= 0; 0 when open), the signed normal gap/clearance `gap` (>0 open, <0 penetrating),
// and the tangential friction traction magnitude `tau`. Emitted to .dat (CSTR block) and
// .frd (CONTACT dataset). Empty on a deck with no *CONTACT PAIR.
struct ContactPoint {
  Index node_id{};
  bool closed{false};  // true = in contact (CLOSED); false = open
  Real p{0.0};         // normal contact pressure (>= 0)
  Real gap{0.0};       // signed normal gap: >0 clearance, <0 penetration
  Real tau{0.0};       // tangential (friction) traction magnitude
};

// Nodal result fields for a static step, aligned with mesh node indices.
struct StaticFields {
  std::vector<Vec3> displacement;  // U
  std::vector<Voigt6> stress;      // S (averaged nodal stress)
  std::vector<Voigt6> strain;      // E (averaged nodal strain, engineering shear)
  std::vector<Vec3> reaction;      // RF (f_int - f_ext)
  // Contact results per slave node (CSTR), populated by the contact driver; empty for a
  // contact-free deck, so the result writers are unchanged without contact. (spec:
  // contact — contact modifiers and output.)
  std::vector<ContactPoint> contact;
};

// Integration-point heat flux q = -k grad(T) at one Gauss point of one element
// (spec: results-output — *EL PRINT HFL). CalculiX prints HFL per element per
// integration point as (qx, qy, qz); `elem_id` is the element's user id and `gp` the
// 1-based integration-point number.
struct HeatFluxPoint {
  Index elem_id{};
  int gp{};        // 1-based integration point index
  Vec3 flux{0, 0, 0};  // q = -k grad(T)
};

// Nodal result fields for a heat-transfer step, aligned with mesh node indices
// (spec: heat-transfer-analysis / results-output). `temperature` is the NT field;
// `flux_reaction` is the concentrated nodal heat-flux reaction RFL at prescribed-
// temperature nodes (Kt*T at the constrained DOFs), analogous to mechanical RF.
// `heat_flux` is the nodal (extrapolated + averaged) heat-flux vector HFL for the
// .frd; `hfl_points` is the raw integration-point HFL for the *EL PRINT HFL .dat.
struct ThermalFields {
  std::vector<Real> temperature;     // NT
  std::vector<Real> flux_reaction;   // RFL
  std::vector<Vec3> heat_flux;       // HFL (nodal, extrapolated) — for .frd
  std::vector<HeatFluxPoint> hfl_points;  // HFL at integration points — for .dat
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
