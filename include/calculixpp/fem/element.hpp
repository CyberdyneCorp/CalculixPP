#pragma once
#include <array>
#include <span>
#include <vector>

#include "calculixpp/core/element.hpp"
#include "calculixpp/core/material.hpp"
#include "calculixpp/core/types.hpp"

// Tetrahedral linear-elastic element kernels.
// Math reference (reimplemented, not copied): CalculiX shape4tet.f / shape10tet.f /
// gauss.f / linel.f / e_c3d.f. Natural coords (xi,et,ze), barycentric, corner-1 at
// origin. Stress/strain order: xx, yy, zz, xy, xz, yz (engineering shear).
namespace cxpp::fem {

inline constexpr int kMaxNodes = 10;
inline constexpr int kVoigt = 6;

struct GaussPoint {
  Real xi, et, ze, w;
};

// Gauss rule used for the stiffness of each element type
// (C3D4: 1 point w=1/6; C3D10: 4 points w=1/24 at (5±√5)/20).
std::span<const GaussPoint> gauss_rule(ElementType type);

// Shape functions and their natural derivatives at a point.
struct Shape {
  int n{};
  std::array<Real, kMaxNodes> N{};
  std::array<std::array<Real, 3>, kMaxNodes> dNdxi{};  // [node][xi,et,ze]
};
Shape shape(ElementType type, Real xi, Real et, Real ze);

// Physical gradients g[k] = dN_k/dX at a point, given node coordinates.
// Fills `grad` (n x 3) and returns det(J). Throws std::runtime_error if det(J) <= 0.
Real physical_gradients(const Shape& s, std::span<const Vec3> coords,
                        std::array<std::array<Real, 3>, kMaxNodes>& grad);

// 6x6 isotropic elastic constitutive matrix (order xx,yy,zz,xy,xz,yz).
using D6 = std::array<std::array<Real, kVoigt>, kVoigt>;
D6 elastic_iso_D(const ElasticIso& mat);

// Element volume via the element's Gauss rule.
Real element_volume(ElementType type, std::span<const Vec3> coords);

// Dense element stiffness (3n x 3n, row-major), DOF order node-major
// [u1,v1,w1, u2,v2,w2, ...]. Ke = Σ_gp Bᵀ D B · det(J) · w.
std::vector<Real> element_stiffness(ElementType type, std::span<const Vec3> coords,
                                    const ElasticIso& mat);

}  // namespace cxpp::fem
