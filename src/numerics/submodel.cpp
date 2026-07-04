#include "calculixpp/numerics/submodel.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>

#include "calculixpp/fem/element.hpp"
#include "calculixpp/numerics/linear_static.hpp"

namespace cxpp::numerics {
namespace {

// Solve the 3x3 system J d = r for d (Cramer's rule). Returns false when J is singular
// (|det| below eps), leaving d untouched — the caller treats that as a failed step.
bool solve3(const std::array<std::array<Real, 3>, 3>& J, const Vec3& r, Vec3& d) {
  const Real c00 = J[1][1] * J[2][2] - J[1][2] * J[2][1];
  const Real c01 = J[1][2] * J[2][0] - J[1][0] * J[2][2];
  const Real c02 = J[1][0] * J[2][1] - J[1][1] * J[2][0];
  const Real det = J[0][0] * c00 + J[0][1] * c01 + J[0][2] * c02;
  if (std::abs(det) < 1e-30) return false;
  const Real inv = 1.0 / det;
  // Column-replacement determinants (Cramer): d_k = det(J with column k <- r) / det.
  const Real dx = r[0] * c00 + r[1] * (J[0][2] * J[2][1] - J[0][1] * J[2][2]) +
                  r[2] * (J[0][1] * J[1][2] - J[0][2] * J[1][1]);
  const Real dy = r[0] * c01 + r[1] * (J[0][0] * J[2][2] - J[0][2] * J[2][0]) +
                  r[2] * (J[0][2] * J[1][0] - J[0][0] * J[1][2]);
  const Real dz = r[0] * c02 + r[1] * (J[0][1] * J[2][0] - J[0][0] * J[2][1]) +
                  r[2] * (J[0][0] * J[1][1] - J[0][1] * J[1][0]);
  d = {dx * inv, dy * inv, dz * inv};
  return true;
}

// Physical image x(ξ) = Σ N_a(ξ) x_a of a natural point in an element with node
// coordinates `coords`.
Vec3 map_forward(const fem::Shape& s, std::span<const Vec3> coords) {
  Vec3 x{0, 0, 0};
  for (int a = 0; a < s.n; ++a) {
    const Real N = s.N[static_cast<std::size_t>(a)];
    for (int k = 0; k < 3; ++k) x[static_cast<std::size_t>(k)] += N * coords[static_cast<std::size_t>(a)][static_cast<std::size_t>(k)];
  }
  return x;
}

// Jacobian J[i][j] = dx_i/dξ_j = Σ_a x_a[i] dN_a/dξ_j at the current natural point.
std::array<std::array<Real, 3>, 3> jacobian(const fem::Shape& s,
                                            std::span<const Vec3> coords) {
  std::array<std::array<Real, 3>, 3> J{};
  for (int a = 0; a < s.n; ++a)
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        J[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] +=
            coords[static_cast<std::size_t>(a)][static_cast<std::size_t>(i)] *
            s.dNdxi[static_cast<std::size_t>(a)][static_cast<std::size_t>(j)];
  return J;
}

// A natural-coordinate seed at the element's reference centroid — the Newton start.
void reference_centroid(ElementType type, Real& xi, Real& et, Real& ze) {
  switch (type) {
    case ElementType::C3D4:
    case ElementType::C3D10:
      xi = et = ze = 0.25;  // tet barycentric centroid
      return;
    case ElementType::C3D6:
    case ElementType::C3D15:
      xi = et = 1.0 / 3.0;  // triangle centroid in (ξ,η)
      ze = 0.0;             // mid-height
      return;
    default:  // hex family: origin of [-1,1]^3
      xi = et = ze = 0.0;
      return;
  }
}

// Is ξ inside the element's reference domain within `tol`? Tets are barycentric
// (ξ,η,ζ ≥ 0, ξ+η+ζ ≤ 1); the hex family is [-1,1]^3; the wedge is the unit triangle in
// (ξ,η) tensored with ζ ∈ [-1,1].
bool inside_domain(ElementType type, Real xi, Real et, Real ze, Real tol) {
  switch (type) {
    case ElementType::C3D4:
    case ElementType::C3D10:
      return xi >= -tol && et >= -tol && ze >= -tol && xi + et + ze <= 1.0 + tol;
    case ElementType::C3D6:
    case ElementType::C3D15:
      return xi >= -tol && et >= -tol && xi + et <= 1.0 + tol && ze >= -1.0 - tol &&
             ze <= 1.0 + tol;
    default:
      return xi >= -1.0 - tol && xi <= 1.0 + tol && et >= -1.0 - tol &&
             et <= 1.0 + tol && ze >= -1.0 - tol && ze <= 1.0 + tol;
  }
}

Real dist(const Vec3& a, const Vec3& b) {
  const Real dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// How far ξ lies OUTSIDE the element's reference domain (0 when inside). For an affine
// element the inverse map reaches any physical point with ~zero residual, so a small
// residual alone does NOT prove containment; the host-element search must also check
// that ξ is at (or barely past) the domain boundary. This measures that overshoot in
// natural-coordinate units.
Real domain_overshoot(ElementType type, Real xi, Real et, Real ze) {
  const auto neg = [](Real v) { return v < 0 ? -v : 0.0; };  // amount below a lower bound
  switch (type) {
    case ElementType::C3D4:
    case ElementType::C3D10:
      return neg(xi) + neg(et) + neg(ze) + neg(1.0 - xi - et - ze);
    case ElementType::C3D6:
    case ElementType::C3D15:
      return neg(xi) + neg(et) + neg(1.0 - xi - et) + neg(1.0 + ze) + neg(1.0 - ze);
    default:  // hex family [-1,1]^3
      return neg(1.0 + xi) + neg(1.0 - xi) + neg(1.0 + et) + neg(1.0 - et) +
             neg(1.0 + ze) + neg(1.0 - ze);
  }
}

}  // namespace

NaturalCoords natural_coords(ElementType type, std::span<const Vec3> coords,
                             const Vec3& X) {
  Real xi = 0, et = 0, ze = 0;
  reference_centroid(type, xi, et, ze);
  Vec3 x{0, 0, 0};
  // Newton: solve x(ξ) = X. The isoparametric map is affine for a linear tet (exact in
  // one step) and mildly nonlinear otherwise; a dozen steps is ample.
  for (int iter = 0; iter < 30; ++iter) {
    const fem::Shape s = fem::shape(type, xi, et, ze);
    x = map_forward(s, coords);
    const Vec3 r{X[0] - x[0], X[1] - x[1], X[2] - x[2]};
    if (dist(x, X) < 1e-10) break;
    Vec3 d{0, 0, 0};
    if (!solve3(jacobian(s, coords), r, d)) break;
    xi += d[0];
    et += d[1];
    ze += d[2];
  }
  NaturalCoords nc;
  nc.xi = xi;
  nc.et = et;
  nc.ze = ze;
  nc.distance = dist(x, X);
  nc.overshoot = domain_overshoot(type, xi, et, ze);
  // A point counts as inside when the converged image sits on X (small residual) AND ξ
  // is within the reference domain. The residual guard rejects a spurious in-domain ξ
  // whose image did not actually reach X (a degenerate/non-convergent map).
  nc.inside = nc.distance < 1e-6 && inside_domain(type, xi, et, ze, 1e-6);
  return nc;
}

namespace {

// Interpolate the global displacement at natural coordinates `nc` within global element
// `e`: Σ_a N_a(ξ) U_a over the element's nodes.
Vec3 interp_at(const GlobalSolution& global, std::size_t e, const NaturalCoords& nc) {
  const std::vector<Index>& conn = global.elem_conn[e];
  const fem::Shape s = fem::shape(global.elem_type[e], nc.xi, nc.et, nc.ze);
  Vec3 u{0, 0, 0};
  for (std::size_t a = 0; a < conn.size(); ++a) {
    const Vec3& ua = global.displacement[static_cast<std::size_t>(conn[a])];
    const Real N = s.N[a];
    for (int k = 0; k < 3; ++k) u[static_cast<std::size_t>(k)] += N * ua[static_cast<std::size_t>(k)];
  }
  return u;
}

}  // namespace

Vec3 interpolate_global_displacement(const GlobalSolution& global, const Vec3& X,
                                     const std::vector<bool>& mask) {
  // Search the (masked) global elements for the host of X: the one whose inverse map
  // places X inside its reference domain. Track the element with the SMALLEST domain
  // overshoot (with X actually on its image) as a fallback for a point sitting a hair
  // outside every element — a cut node landing exactly on a global element face, nudged
  // out by round-off. A genuinely external point has a large overshoot AND is rejected.
  Real best_overshoot = std::numeric_limits<Real>::max();
  std::size_t best_elem = global.num_elements();
  NaturalCoords best_nc;
  for (std::size_t e = 0; e < global.num_elements(); ++e) {
    if (!mask.empty() && !mask[e]) continue;
    const std::vector<Index>& conn = global.elem_conn[e];
    std::vector<Vec3> coords(conn.size());
    for (std::size_t a = 0; a < conn.size(); ++a)
      coords[a] = global.coords[static_cast<std::size_t>(conn[a])];
    const NaturalCoords nc = natural_coords(global.elem_type[e], coords, X);
    if (nc.inside) return interp_at(global, e, nc);
    // Only elements whose image actually reaches X (small residual) are fallback
    // candidates; among them keep the one X is least outside.
    if (nc.distance < 1e-6 && nc.overshoot < best_overshoot) {
      best_overshoot = nc.overshoot;
      best_elem = e;
      best_nc = nc;
    }
  }
  // Accept the closest element only when X is barely outside it (round-off on a shared
  // face). Beyond that natural-coordinate tolerance, no host element contains X.
  if (best_elem < global.num_elements() && best_overshoot < 1e-6)
    return interp_at(global, best_elem, best_nc);
  throw std::runtime_error(
      "submodel: no host global element found for boundary node at (" +
      std::to_string(X[0]) + ", " + std::to_string(X[1]) + ", " +
      std::to_string(X[2]) + "); it falls outside the global element set");
}

StaticFields solve_submodel(const Model& model, const GlobalSolution& global,
                            std::optional<compute::SolverKind> forced) {
  // The global element set of the first *SUBMODEL card bounds the host search (empty ->
  // every global element). Build the per-global-element mask once.
  std::vector<bool> mask;  // empty -> search all global elements
  if (!model.submodels.empty() && !model.submodels.front().global_elset.empty()) {
    // The global elset is expressed in GLOBAL element ids; a GlobalSolution carries no
    // ids, so the caller is expected to pass a GlobalSolution already restricted to the
    // set. We keep the mask empty here (search all supplied global elements) — the
    // supplied global mesh IS the searched set. (Documented in tasks: global-elset
    // resolution against a full global deck is part of the .frd follow-on.)
  }

  // Fill the driven SPCs with the interpolated global displacement at each node.
  Model driven = model;
  for (Spc& spc : driven.spcs) {
    if (!spc.driven) continue;
    const Index ni = driven.mesh.node_index(spc.node_id);
    if (ni < 0)
      throw std::runtime_error("submodel: driven boundary node " +
                               std::to_string(spc.node_id) + " not in the submodel mesh");
    const Vec3& X = driven.mesh.nodes()[static_cast<std::size_t>(ni)].x;
    const Vec3 u = interpolate_global_displacement(global, X, mask);
    spc.value = u[static_cast<std::size_t>(spc.comp - 1)];
    spc.driven = false;  // now a literal prescribed value for the ordinary solve
  }
  return solve_linear_static(driven, forced);
}

}  // namespace cxpp::numerics
