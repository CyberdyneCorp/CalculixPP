#pragma once
#include <array>
#include <cstdint>

namespace cxpp {

// Central index type (spec: phase-1 design decision — 32-bit signed `cpp_index_t`,
// covers mobile-scale models < 2.1B DOF; widen here if ever needed).
using Index = std::int32_t;
using Real = double;
using Vec3 = std::array<Real, 3>;

inline constexpr int kDofsPerNode = 3;  // structural: u, v, w

}  // namespace cxpp
