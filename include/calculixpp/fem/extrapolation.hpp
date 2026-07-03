#pragma once
#include <array>
#include <cstddef>
#include <vector>

#include "calculixpp/core/element.hpp"
#include "calculixpp/core/types.hpp"

// Integration-point -> nodal extrapolation, shared by stress (Voigt6) recovery and
// heat-flux (Vec3) recovery so both use the identical scheme (spec: results-output —
// reimplemented from CalculiX extrapolate.f, not copied). C3D4 is constant; C3D10
// uses the a4 corner extrapolation matrix + nonei10 midside averaging; the hex/wedge
// families use the Gauss-point mean at every node (an order-0 recovery — the nodal
// value is not part of the displacement-accuracy gate, and CalculiX itself reports HFL
// at the integration points, which we emit verbatim to the .dat).
namespace cxpp::fem {

// C3D10 corner extrapolation matrix a4 (inverse of the linear-tet shape functions at
// the 4 Gauss points) and midside adjacency (nonei10), from extrapolate.f.
inline constexpr Real kA4diag = 1.927050983124842;
inline constexpr Real kA4off = -0.309016994374947;
inline constexpr int kMid[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};

// Extrapolate `NC`-component integration-point values `gp` to the `n` element nodes.
// `Comp` is a fixed-size array-like type (std::array<Real,NC>) indexed by [0..NC).
template <std::size_t NC>
void extrapolate_gp_to_nodes(ElementType type, int n,
                             const std::vector<std::array<Real, NC>>& gp,
                             std::vector<std::array<Real, NC>>& nodal) {
  if (type == ElementType::C3D4) {
    for (int i = 0; i < n; ++i) nodal[static_cast<std::size_t>(i)] = gp[0];
    return;
  }
  if (type != ElementType::C3D10) {
    std::array<Real, NC> mean{};
    for (const auto& v : gp)
      for (std::size_t c = 0; c < NC; ++c) mean[c] += v[c];
    const Real inv = gp.empty() ? 0.0 : 1.0 / static_cast<Real>(gp.size());
    for (std::size_t c = 0; c < NC; ++c) mean[c] *= inv;
    for (int i = 0; i < n; ++i) nodal[static_cast<std::size_t>(i)] = mean;
    return;
  }
  for (int c = 0; c < 4; ++c)
    for (std::size_t comp = 0; comp < NC; ++comp) {
      Real v = 0.0;
      for (int l = 0; l < 4; ++l)
        v += (c == l ? kA4diag : kA4off) * gp[static_cast<std::size_t>(l)][comp];
      nodal[static_cast<std::size_t>(c)][comp] = v;
    }
  for (int mnode = 0; mnode < 6; ++mnode)
    for (std::size_t comp = 0; comp < NC; ++comp)
      nodal[static_cast<std::size_t>(4 + mnode)][comp] =
          0.5 * (nodal[static_cast<std::size_t>(kMid[mnode][0])][comp] +
                 nodal[static_cast<std::size_t>(kMid[mnode][1])][comp]);
}

}  // namespace cxpp::fem
