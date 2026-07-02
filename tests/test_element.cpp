// Element-math unit tests: shape functions, physical gradients (consistency),
// element volume, and stiffness properties (symmetry, rigid-body null space).
#include <array>
#include <vector>

#include "calculixpp/fem/element.hpp"
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

}  // namespace

int main() {
  test_partition_of_unity();
  test_gradients_and_volume();
  test_stiffness_properties();
  if (cxtest::g_failures == 0) std::printf("test_element: OK\n");
  CX_MAIN_RETURN();
}
