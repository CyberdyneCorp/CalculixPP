// Element-math unit tests: shape functions, physical gradients (consistency),
// element volume, and stiffness properties (symmetry, rigid-body null space).
#include <array>
#include <vector>

#include "calculixpp/fem/element.hpp"
#include "calculixpp/fem/material_model.hpp"
#include "check.hpp"

using namespace cxpp;
using namespace cxpp::fem;

namespace {

// Build a straight-edged C3D10 (edge nodes at corner midpoints, per nonei10 order).
std::vector<Vec3> c3d10_from_corners(const std::array<Vec3, 4>& c) {
  auto mid = [](const Vec3& a, const Vec3& b) {
    return Vec3{(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2};
  };
  return {c[0], c[1], c[2], c[3],
          mid(c[0], c[1]),  // N5  edge 1-2
          mid(c[1], c[2]),  // N6  edge 2-3
          mid(c[2], c[0]),  // N7  edge 3-1
          mid(c[0], c[3]),  // N8  edge 1-4
          mid(c[1], c[3]),  // N9  edge 2-4
          mid(c[2], c[3])}; // N10 edge 3-4
}

void test_partition_of_unity() {
  for (ElementType t : {ElementType::C3D4, ElementType::C3D10}) {
    const Shape s = shape(t, 0.2, 0.3, 0.15);
    Real sumN = 0.0;
    std::array<Real, 3> sumd{0, 0, 0};
    for (int k = 0; k < s.n; ++k) {
      sumN += s.N[static_cast<std::size_t>(k)];
      for (int d = 0; d < 3; ++d) sumd[static_cast<std::size_t>(d)] += s.dNdxi[static_cast<std::size_t>(k)][static_cast<std::size_t>(d)];
    }
    CX_NEAR(sumN, 1.0, 1e-12);
    for (int d = 0; d < 3; ++d) CX_NEAR(sumd[static_cast<std::size_t>(d)], 0.0, 1e-12);
  }
}

// Σ_k g_k = 0 and Σ_k g_k ⊗ x_k = I — rigorous check of the Jacobian/transpose.
void check_gradient_consistency(ElementType t, const std::vector<Vec3>& coords) {
  std::array<std::array<Real, 3>, kMaxNodes> g{};
  for (const GaussPoint& gp : gauss_rule(t)) {
    const Shape s = shape(t, gp.xi, gp.et, gp.ze);
    const Real det = physical_gradients(s, coords, g);
    CX_CHECK(det > 0.0);
    std::array<Real, 3> sg{0, 0, 0};
    std::array<std::array<Real, 3>, 3> repro{};
    for (int k = 0; k < s.n; ++k)
      for (int i = 0; i < 3; ++i) {
        sg[static_cast<std::size_t>(i)] += g[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)];
        for (int j = 0; j < 3; ++j)
          repro[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] +=
              g[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)] * coords[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)];
      }
    for (int i = 0; i < 3; ++i) {
      CX_NEAR(sg[static_cast<std::size_t>(i)], 0.0, 1e-10);
      for (int j = 0; j < 3; ++j) CX_NEAR(repro[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)], i == j ? 1.0 : 0.0, 1e-10);
    }
  }
}

void test_gradients_and_volume() {
  const std::vector<Vec3> tet4{{0, 0, 0}, {2, 0, 0}, {1, 3, 0}, {0.5, 0.7, 4}};
  check_gradient_consistency(ElementType::C3D4, tet4);

  const std::array<Vec3, 4> unit{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
  const std::vector<Vec3> tet10 = c3d10_from_corners(unit);
  check_gradient_consistency(ElementType::C3D10, tet10);

  // Unit tet volume = 1/6 for both element types.
  const std::vector<Vec3> unit4{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  CX_NEAR(element_volume(ElementType::C3D4, unit4), 1.0 / 6.0, 1e-12);
  CX_NEAR(element_volume(ElementType::C3D10, tet10), 1.0 / 6.0, 1e-12);
}

void test_stiffness_properties() {
  const ElasticIso mat{210000.0, 0.3};
  const std::vector<Vec3> coords{{0, 0, 0}, {2, 0, 0}, {1, 3, 0}, {0.5, 0.7, 4}};
  const std::vector<Real> Ke = element_stiffness(ElementType::C3D4, coords, mat);
  const int ndof = 12;

  // Symmetry.
  for (int i = 0; i < ndof; ++i)
    for (int j = 0; j < ndof; ++j)
      CX_NEAR(Ke[static_cast<std::size_t>(i * ndof + j)], Ke[static_cast<std::size_t>(j * ndof + i)], 1e-6);

  // Rigid-body translations lie in the null space: Ke * t ≈ 0.
  for (int dir = 0; dir < 3; ++dir) {
    std::array<Real, 12> t{};
    for (int nd = 0; nd < 4; ++nd) t[static_cast<std::size_t>(nd * 3 + dir)] = 1.0;
    for (int i = 0; i < ndof; ++i) {
      Real r = 0.0;
      for (int j = 0; j < ndof; ++j) r += Ke[static_cast<std::size_t>(i * ndof + j)] * t[static_cast<std::size_t>(j)];
      CX_NEAR(r, 0.0, 1e-6);
    }
  }
}

// The material-point kernel (element_tangent_force with ElasticIsoMaterial) must
// reproduce the closed-form element_stiffness tangent, and give fe = Ke*ue for an
// arbitrary displacement (small-strain internal force consistency).
void test_material_point_kernel(ElementType type, const std::vector<Vec3>& coords) {
  const ElasticIso mat{210000.0, 0.3};
  const int n = nodes_per_element(type);
  const int ndof = n * 3;
  const std::vector<Real> Ke_ref = element_stiffness(type, coords, mat);

  ElasticIsoMaterial material(mat);
  std::vector<MaterialState> state(gauss_rule(type).size());
  // Tangent at zero displacement must match element_stiffness to ~1e-12.
  const std::vector<Vec3> zero(static_cast<std::size_t>(n), Vec3{0, 0, 0});
  const ElementResponse r0 = element_tangent_force(type, coords, zero, material, state);
  for (int i = 0; i < ndof * ndof; ++i)
    CX_NEAR(r0.Ke[static_cast<std::size_t>(i)], Ke_ref[static_cast<std::size_t>(i)],
            1e-12 * (1.0 + std::fabs(Ke_ref[static_cast<std::size_t>(i)])));
  for (int i = 0; i < ndof; ++i) CX_NEAR(r0.fe[static_cast<std::size_t>(i)], 0.0, 1e-9);

  // For a nonzero displacement, fe == Ke * ue (linear elasticity), and Ke unchanged.
  std::vector<Vec3> ue(static_cast<std::size_t>(n));
  for (int k = 0; k < n; ++k)
    ue[static_cast<std::size_t>(k)] = {0.001 * (k + 1), -0.002 * k, 0.0007 * (k + 2)};
  std::vector<MaterialState> state2(gauss_rule(type).size());
  const ElementResponse r1 = element_tangent_force(type, coords, ue, material, state2);
  for (int a = 0; a < ndof; ++a) {
    Real ku = 0.0;
    for (int b = 0; b < ndof; ++b)
      ku += Ke_ref[static_cast<std::size_t>(a * ndof + b)] *
            ue[static_cast<std::size_t>(b / 3)][static_cast<std::size_t>(b % 3)];
    CX_NEAR(r1.fe[static_cast<std::size_t>(a)], ku, 1e-7 * (1.0 + std::fabs(ku)));
  }
}

void test_material_point() {
  const std::vector<Vec3> tet4{{0, 0, 0}, {2, 0, 0}, {1, 3, 0}, {0.5, 0.7, 4}};
  test_material_point_kernel(ElementType::C3D4, tet4);
  const std::array<Vec3, 4> unit{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
  test_material_point_kernel(ElementType::C3D10, c3d10_from_corners(unit));
}

// ---- Hex/wedge element families (workstream 3.2/3.3) --------------------------

// Corner coordinates of a unit hexahedron in CalculiX C3D8 node order.
const std::array<Vec3, 8> kHexCorners{{
    {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
    {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1},
}};

// Build a C3D8 element from 8 corners.
std::vector<Vec3> c3d8_coords(const std::array<Vec3, 8>& c) {
  return {c.begin(), c.end()};
}

// Build a straight-edged C3D20 from 8 corners (edge midsides at corner midpoints),
// C3D20 edge ordering: 9-12 bottom, 13-16 top, 17-20 verticals.
std::vector<Vec3> c3d20_coords(const std::array<Vec3, 8>& c) {
  auto mid = [](const Vec3& a, const Vec3& b) {
    return Vec3{(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2};
  };
  const int e[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
                        {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
  std::vector<Vec3> out(c.begin(), c.end());
  for (auto& pr : e) out.push_back(mid(c[static_cast<std::size_t>(pr[0])],
                                       c[static_cast<std::size_t>(pr[1])]));
  return out;
}

// Wedge corners (C3D6): triangle base z=0 then z=1.
const std::array<Vec3, 6> kWedgeCorners{{
    {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1},
}};
std::vector<Vec3> c3d6_coords(const std::array<Vec3, 6>& c) {
  return {c.begin(), c.end()};
}
// C3D15: 6 corners then 9 midsides (edges 1-2,2-3,3-1 bottom; 4-5,5-6,6-4 top;
// 1-4,2-5,3-6 verticals) per shape15w.f ordering 7..15.
std::vector<Vec3> c3d15_coords(const std::array<Vec3, 6>& c) {
  auto mid = [](const Vec3& a, const Vec3& b) {
    return Vec3{(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2};
  };
  const int e[9][2] = {{0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5},
                       {5, 3}, {0, 3}, {1, 4}, {2, 5}};
  std::vector<Vec3> out(c.begin(), c.end());
  for (auto& pr : e) out.push_back(mid(c[static_cast<std::size_t>(pr[0])],
                                       c[static_cast<std::size_t>(pr[1])]));
  return out;
}

// Partition of unity + zero derivative sum at a few sample points.
void check_partition_of_unity(ElementType t) {
  for (const GaussPoint& gp : gauss_rule(t)) {
    const Shape s = shape(t, gp.xi, gp.et, gp.ze);
    Real sumN = 0.0;
    std::array<Real, 3> sumd{0, 0, 0};
    for (int k = 0; k < s.n; ++k) {
      sumN += s.N[static_cast<std::size_t>(k)];
      for (int d = 0; d < 3; ++d)
        sumd[static_cast<std::size_t>(d)] += s.dNdxi[static_cast<std::size_t>(k)][static_cast<std::size_t>(d)];
    }
    CX_NEAR(sumN, 1.0, 1e-12);
    for (int d = 0; d < 3; ++d) CX_NEAR(sumd[static_cast<std::size_t>(d)], 0.0, 1e-12);
  }
}

// Gauss-rule volume must equal the analytic element volume.
void check_volume(ElementType t, const std::vector<Vec3>& coords, Real expected) {
  CX_NEAR(element_volume(t, coords), expected, 1e-12);
}

// Constant-strain patch test: impose a linear field u_k = eps . x_k on the nodes and
// verify the recovered strain equals eps at every Gauss point (a superconvergent
// isoparametric-element property, and the key correctness check for a new topology).
void check_patch(ElementType t, const std::vector<Vec3>& coords) {
  const int n = nodes_per_element(t);
  // A general symmetric small strain (engineering shear in xy,xz,yz).
  const Real exx = 1e-3, eyy = -4e-4, ezz = 2e-4, gxy = 6e-4, gxz = -3e-4, gyz = 5e-4;
  std::vector<Vec3> ue(static_cast<std::size_t>(n));
  for (int k = 0; k < n; ++k) {
    const Real x = coords[static_cast<std::size_t>(k)][0];
    const Real y = coords[static_cast<std::size_t>(k)][1];
    const Real z = coords[static_cast<std::size_t>(k)][2];
    ue[static_cast<std::size_t>(k)] = {exx * x + 0.5 * gxy * y + 0.5 * gxz * z,
                                       0.5 * gxy * x + eyy * y + 0.5 * gyz * z,
                                       0.5 * gxz * x + 0.5 * gyz * y + ezz * z};
  }
  std::array<std::array<Real, 3>, kMaxNodes> g{};
  for (const GaussPoint& gp : gauss_rule(t)) {
    const Shape s = shape(t, gp.xi, gp.et, gp.ze);
    const Real det = physical_gradients(s, coords, g);
    CX_CHECK(det > 0.0);
    const Voigt6 e = strain_from_gradients(g, n, ue);
    CX_NEAR(e[0], exx, 1e-11);
    CX_NEAR(e[1], eyy, 1e-11);
    CX_NEAR(e[2], ezz, 1e-11);
    CX_NEAR(e[3], gxy, 1e-11);
    CX_NEAR(e[4], gxz, 1e-11);
    CX_NEAR(e[5], gyz, 1e-11);
  }
}

// Stiffness symmetry and the six rigid-body modes (3 translations + 3 rotations)
// must lie in the null space of a full-integration element.
void check_rigid_body(ElementType t, const std::vector<Vec3>& coords) {
  const ElasticIso mat{210000.0, 0.3};
  const int n = nodes_per_element(t);
  const int ndof = n * 3;
  const std::vector<Real> Ke = element_stiffness(t, coords, mat);
  for (int i = 0; i < ndof; ++i)
    for (int j = 0; j < ndof; ++j)
      CX_NEAR(Ke[static_cast<std::size_t>(i * ndof + j)],
              Ke[static_cast<std::size_t>(j * ndof + i)], 1e-6);
  // Rigid-body fields: 3 translations and 3 infinitesimal rotations.
  auto rbm = [&](int m) {
    std::vector<Real> v(static_cast<std::size_t>(ndof), 0.0);
    for (int k = 0; k < n; ++k) {
      const Real x = coords[static_cast<std::size_t>(k)][0];
      const Real y = coords[static_cast<std::size_t>(k)][1];
      const Real z = coords[static_cast<std::size_t>(k)][2];
      Vec3 d{0, 0, 0};
      if (m == 0) d = {1, 0, 0};
      else if (m == 1) d = {0, 1, 0};
      else if (m == 2) d = {0, 0, 1};
      else if (m == 3) d = {-y, x, 0};   // rot z
      else if (m == 4) d = {-z, 0, x};   // rot y
      else d = {0, -z, y};               // rot x
      for (int c = 0; c < 3; ++c) v[static_cast<std::size_t>(k * 3 + c)] = d[static_cast<std::size_t>(c)];
    }
    return v;
  };
  for (int m = 0; m < 6; ++m) {
    const std::vector<Real> v = rbm(m);
    Real norm = 0.0;
    for (int i = 0; i < ndof; ++i) {
      Real r = 0.0;
      for (int j = 0; j < ndof; ++j)
        r += Ke[static_cast<std::size_t>(i * ndof + j)] * v[static_cast<std::size_t>(j)];
      norm += r * r;
    }
    CX_NEAR(std::sqrt(norm), 0.0, 1e-3);
  }
}

// The C3D8R hourglass control must give a strictly positive strain energy to a
// pure hourglass displacement (mode Γ1 applied in x), which the un-stabilized
// single-point element would leave at zero energy.
void check_hourglass_suppressed(const std::vector<Vec3>& h8) {
  const ElasticIso mat{210000.0, 0.3};
  const std::vector<Real> Ke = element_stiffness(ElementType::C3D8R, h8, mat);
  // Hourglass base vector Γ1 (per corner), applied to the x-displacement.
  const Real gamma1[8] = {1, 1, -1, -1, -1, -1, 1, 1};
  std::array<Real, 24> u{};
  for (int k = 0; k < 8; ++k) u[static_cast<std::size_t>(k * 3)] = gamma1[k];
  Real energy = 0.0;
  for (int i = 0; i < 24; ++i)
    for (int j = 0; j < 24; ++j)
      energy += u[static_cast<std::size_t>(i)] *
                Ke[static_cast<std::size_t>(i * 24 + j)] * u[static_cast<std::size_t>(j)];
  CX_CHECK(energy > 1.0);  // strictly positive -> mode is suppressed
}

void test_hex_wedge() {
  const std::vector<Vec3> h8 = c3d8_coords(kHexCorners);
  const std::vector<Vec3> h20 = c3d20_coords(kHexCorners);
  const std::vector<Vec3> w6 = c3d6_coords(kWedgeCorners);
  const std::vector<Vec3> w15 = c3d15_coords(kWedgeCorners);

  for (ElementType t : {ElementType::C3D8, ElementType::C3D8R, ElementType::C3D20,
                        ElementType::C3D20R, ElementType::C3D6, ElementType::C3D15})
    check_partition_of_unity(t);

  check_volume(ElementType::C3D8, h8, 1.0);
  check_volume(ElementType::C3D20, h20, 1.0);
  check_volume(ElementType::C3D6, w6, 0.5);
  check_volume(ElementType::C3D15, w15, 0.5);

  // Patch tests: constant strain recovered exactly at every Gauss point (also for
  // the reduced elements — the centroid/2x2x2 rules integrate a constant field
  // exactly, and hourglass control is orthogonal to the linear field).
  check_patch(ElementType::C3D8, h8);
  check_patch(ElementType::C3D8R, h8);
  check_patch(ElementType::C3D20, h20);
  check_patch(ElementType::C3D20R, h20);
  check_patch(ElementType::C3D6, w6);
  check_patch(ElementType::C3D15, w15);

  // Rigid-body null space. Includes the reduced elements: C3D8R must keep exactly
  // the 6 rigid-body modes as zero-energy (hourglass control removes the 12
  // spurious modes but must NOT stiffen the rigid-body ones).
  check_rigid_body(ElementType::C3D8, h8);
  check_rigid_body(ElementType::C3D8R, h8);
  check_rigid_body(ElementType::C3D20, h20);
  check_rigid_body(ElementType::C3D20R, h20);
  check_rigid_body(ElementType::C3D6, w6);
  check_rigid_body(ElementType::C3D15, w15);

  // Hourglass control must suppress the spurious modes: an hourglass displacement
  // (mode Γ1 in x) produces zero strain at the centroid, so the un-stabilized
  // reduced stiffness gives ~0 energy; with hourglass control the energy is > 0.
  check_hourglass_suppressed(h8);

  // Material-point kernel reproduces element_stiffness (incl. hourglass) and gives
  // fe = Ke ue for every family.
  test_material_point_kernel(ElementType::C3D8, h8);
  test_material_point_kernel(ElementType::C3D8R, h8);
  test_material_point_kernel(ElementType::C3D20, h20);
  test_material_point_kernel(ElementType::C3D20R, h20);
  test_material_point_kernel(ElementType::C3D6, w6);
  test_material_point_kernel(ElementType::C3D15, w15);
}

}  // namespace

int main() {
  test_partition_of_unity();
  test_gradients_and_volume();
  test_stiffness_properties();
  test_material_point();
  test_hex_wedge();
  if (cxtest::g_failures == 0) std::printf("test_element: OK\n");
  CX_MAIN_RETURN();
}
