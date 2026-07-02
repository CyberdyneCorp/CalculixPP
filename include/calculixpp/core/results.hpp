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
  std::vector<Vec3> reaction;      // RF (f_int - f_ext)
};

}  // namespace cxpp
