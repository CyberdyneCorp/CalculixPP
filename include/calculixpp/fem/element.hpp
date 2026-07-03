#pragma once
#include <array>
#include <span>
#include <vector>

#include "calculixpp/core/element.hpp"
#include "calculixpp/core/material.hpp"
#include "calculixpp/core/results.hpp"
#include "calculixpp/core/types.hpp"
#include "calculixpp/fem/material_model.hpp"

// Tetrahedral linear-elastic element kernels.
// Math reference (reimplemented, not copied): CalculiX shape4tet.f / shape10tet.f /
// gauss.f / linel.f / e_c3d.f. Natural coords (xi,et,ze), barycentric, corner-1 at
// origin. Stress/strain order: xx, yy, zz, xy, xz, yz (engineering shear).
namespace cxpp::fem {

inline constexpr int kMaxNodes = 20;
inline constexpr int kVoigt = kVoigt6;  // Voigt dimension (see material_model.hpp)

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
// (D6 is defined in material_model.hpp.)
D6 elastic_iso_D(const ElasticIso& mat);

// Element volume via the element's Gauss rule.
Real element_volume(ElementType type, std::span<const Vec3> coords);

// Strain-displacement entry B[row][3*k+dir] for node k, direction dir (0..2), given
// the physical gradients g (n x 3) of the shape functions. Voigt row order
// xx,yy,zz,xy,xz,yz with engineering shear. Topology-independent: the topology only
// enters through g. (Shared by the stiffness and material-point kernels.)
Real b_entry(int row, const std::array<Real, 3>& gk, int dir);

// Small-strain Voigt6 at physical gradients g for element displacements ue (n nodes).
// strain = B ue with engineering shear (order xx,yy,zz,xy,xz,yz).
Voigt6 strain_from_gradients(const std::array<std::array<Real, 3>, kMaxNodes>& g,
                             int n, std::span<const Vec3> ue);

// Dense element stiffness (3n x 3n, row-major), DOF order node-major
// [u1,v1,w1, u2,v2,w2, ...]. Ke = Σ_gp Bᵀ D B · det(J) · w.
std::vector<Real> element_stiffness(ElementType type, std::span<const Vec3> coords,
                                    const ElasticIso& mat);

// Dense scalar conduction matrix (n x n, row-major) for a temperature field with
// one DOF per node (spec: heat-transfer-analysis — steady conduction). For an
// isotropic conductivity k:
//   Kt_e[a][b] = Σ_gp (g_a · g_b) · k · det(J) · w
// where g_k are the physical shape-function gradients. This reuses the exact same
// shape functions, Gauss rules, and physical_gradients as the mechanical stiffness,
// so every solid element family is supported. (Note the topology-agnostic gradient
// dot product — no B-matrix, since the unknown is scalar.)
std::vector<Real> element_conduction(ElementType type, std::span<const Vec3> coords,
                                     Real k);

// Dense scalar capacitance matrix (n x n, row-major) for the transient thermal term:
//   C_e[a][b] = Σ_gp rho*c · N_a N_b · det(J) · w.
// Zero rho*c yields a zero matrix (steady-state). (spec: heat-transfer-analysis —
// transient conduction; provided now so the transient driver can reuse it.)
std::vector<Real> element_capacitance(ElementType type, std::span<const Vec3> coords,
                                      Real rho_c);

// Thermal-strain equivalent nodal load (spec: heat-transfer — thermal expansion).
// For a thermal strain eps_th = alpha (T - Tref) on the normal components (Voigt
// xx,yy,zz), the mechanical residual gains f_th = ∫ Bᵀ D eps_th dV, so heating a
// constrained body induces the correct thermal stress. `te` is the per-node
// temperature CHANGE (T - Tref) of the element (size n); it is interpolated to each
// Gauss point through the shape functions and multiplied by `alpha`. Returns the 3n
// force vector (DOF order [u1,v1,w1,...]); add it to the mechanical rhs. For zero
// temperature change (or alpha 0) the vector is zero, so the mechanical path is
// unchanged. (Reuses the same Gauss rule / physical_gradients as element_stiffness.)
std::vector<Real> element_thermal_load(ElementType type,
                                       std::span<const Vec3> coords, const D6& D,
                                       Real alpha, std::span<const Real> te);

// Material-point element assembly. For each Gauss point it forms the small-strain
// B ue, evaluates `material` (advancing per-point `state`), and accumulates the
// consistent tangent and internal force:
//   Ke = Σ_gp Bᵀ D_t B · det(J) · w      (3n x 3n, row-major, symmetric layout)
//   fe = Σ_gp Bᵀ σ    · det(J) · w        (size 3n)
// `state` has one entry per Gauss point (gauss_rule(type).size()); pass a vector of
// default MaterialState for the virgin material. For a linear-elastic material and
// ue = 0 this reproduces element_stiffness exactly, and fe = Ke·ue in general.
// (spec: nonlinear-solution-control / material-models — material-point assembly.)
struct ElementResponse {
  std::vector<Real> Ke;  // (3n)^2 row-major
  std::vector<Real> fe;  // 3n
};
ElementResponse element_tangent_force(ElementType type,
                                      std::span<const Vec3> coords,
                                      std::span<const Vec3> ue,
                                      const MaterialModel& material,
                                      std::span<MaterialState> state);

}  // namespace cxpp::fem
