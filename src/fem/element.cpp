#include "calculixpp/fem/element.hpp"

#include <array>
#include <cmath>
#include <stdexcept>

namespace cxpp::fem {
namespace {

// (5 ± √5) / 20 abscissae for the 4-point tetrahedral rule.
constexpr Real kA = 0.138196601125010515;  // (5 - √5)/20
constexpr Real kB = 0.585410196624968454;  // (5 + √5)/20

constexpr std::array<GaussPoint, 1> kGaussC3D4{{{0.25, 0.25, 0.25, 1.0 / 6.0}}};
constexpr std::array<GaussPoint, 4> kGaussC3D10{{
    {kA, kA, kA, 1.0 / 24.0},
    {kB, kA, kA, 1.0 / 24.0},
    {kA, kB, kA, 1.0 / 24.0},
    {kA, kA, kB, 1.0 / 24.0},
}};

// Inverse-transpose of a 3x3 matrix; returns det. J[i][j] = dX_i/dxi_j.
Real invert3_transpose(const std::array<std::array<Real, 3>, 3>& J,
                       std::array<std::array<Real, 3>, 3>& invT) {
  const Real c00 = J[1][1] * J[2][2] - J[1][2] * J[2][1];
  const Real c01 = J[1][2] * J[2][0] - J[1][0] * J[2][2];
  const Real c02 = J[1][0] * J[2][1] - J[1][1] * J[2][0];
  const Real det = J[0][0] * c00 + J[0][1] * c01 + J[0][2] * c02;
  if (det <= 0.0) {
    throw std::runtime_error("non-positive element Jacobian (bad connectivity?)");
  }
  const Real inv = 1.0 / det;
  // J^{-1} rows from cofactors; we want J^{-T} so that g_k = J^{-T} dN_k/dxi.
  const std::array<std::array<Real, 3>, 3> Jinv{{
      {c00 * inv, (J[0][2] * J[2][1] - J[0][1] * J[2][2]) * inv,
       (J[0][1] * J[1][2] - J[0][2] * J[1][1]) * inv},
      {c01 * inv, (J[0][0] * J[2][2] - J[0][2] * J[2][0]) * inv,
       (J[0][2] * J[1][0] - J[0][0] * J[1][2]) * inv},
      {c02 * inv, (J[0][1] * J[2][0] - J[0][0] * J[2][1]) * inv,
       (J[0][0] * J[1][1] - J[0][1] * J[1][0]) * inv},
  }};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) invT[i][j] = Jinv[j][i];
  return det;
}

}  // namespace

std::span<const GaussPoint> gauss_rule(ElementType type) {
  if (type == ElementType::C3D4) return kGaussC3D4;
  return kGaussC3D10;
}

Shape shape(ElementType type, Real xi, Real et, Real ze) {
  Shape s;
  if (type == ElementType::C3D4) {
    s.n = 4;
    s.N = {1.0 - xi - et - ze, xi, et, ze};
    s.dNdxi[0] = {-1, -1, -1};
    s.dNdxi[1] = {1, 0, 0};
    s.dNdxi[2] = {0, 1, 0};
    s.dNdxi[3] = {0, 0, 1};
    return s;
  }
  // C3D10 (shape10tet.f). a = 1 - xi - et - ze.
  const Real a = 1.0 - xi - et - ze;
  s.n = 10;
  s.N = {(2 * a - 1) * a, xi * (2 * xi - 1), et * (2 * et - 1), ze * (2 * ze - 1),
         4 * xi * a,      4 * xi * et,       4 * et * a,        4 * ze * a,
         4 * xi * ze,     4 * et * ze};
  s.dNdxi[0] = {1 - 4 * a, 1 - 4 * a, 1 - 4 * a};
  s.dNdxi[1] = {4 * xi - 1, 0, 0};
  s.dNdxi[2] = {0, 4 * et - 1, 0};
  s.dNdxi[3] = {0, 0, 4 * ze - 1};
  s.dNdxi[4] = {4 * (a - xi), -4 * xi, -4 * xi};
  s.dNdxi[5] = {4 * et, 4 * xi, 0};
  s.dNdxi[6] = {-4 * et, 4 * (a - et), -4 * et};
  s.dNdxi[7] = {-4 * ze, -4 * ze, 4 * (a - ze)};
  s.dNdxi[8] = {4 * ze, 0, 4 * xi};
  s.dNdxi[9] = {0, 4 * ze, 4 * et};
  return s;
}

Real physical_gradients(const Shape& s, std::span<const Vec3> coords,
                        std::array<std::array<Real, 3>, kMaxNodes>& grad) {
  // J[i][j] = sum_k coords[k][i] * dN_k/dxi_j
  std::array<std::array<Real, 3>, 3> J{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      Real v = 0.0;
      for (int k = 0; k < s.n; ++k) v += coords[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)] * s.dNdxi[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)];
      J[i][j] = v;
    }
  std::array<std::array<Real, 3>, 3> invT{};
  const Real det = invert3_transpose(J, invT);
  for (int k = 0; k < s.n; ++k)
    for (int i = 0; i < 3; ++i) {
      Real v = 0.0;
      for (int j = 0; j < 3; ++j) v += invT[i][j] * s.dNdxi[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)];
      grad[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)] = v;
    }
  return det;
}

D6 elastic_iso_D(const ElasticIso& mat) {
  const Real E = mat.E;
  const Real nu = mat.nu;
  const Real um2 = E / (1.0 + nu);           // 2 mu
  const Real lam = nu * um2 / (1.0 - 2.0 * nu);
  const Real mu = um2 / 2.0;
  const Real lam2mu = lam + um2;
  D6 D{};
  D[0][0] = D[1][1] = D[2][2] = lam2mu;
  D[0][1] = D[0][2] = D[1][0] = D[1][2] = D[2][0] = D[2][1] = lam;
  D[3][3] = D[4][4] = D[5][5] = mu;
  return D;
}

Real element_volume(ElementType type, std::span<const Vec3> coords) {
  Real vol = 0.0;
  std::array<std::array<Real, 3>, kMaxNodes> grad{};
  for (const GaussPoint& gp : gauss_rule(type)) {
    const Shape s = shape(type, gp.xi, gp.et, gp.ze);
    const Real det = physical_gradients(s, coords, grad);
    vol += det * gp.w;
  }
  return vol;
}

std::vector<Real> element_stiffness(ElementType type, std::span<const Vec3> coords,
                                    const ElasticIso& mat) {
  const int n = nodes_per_element(type);
  const int ndof = n * kDofsPerNode;
  const D6 D = elastic_iso_D(mat);
  std::vector<Real> Ke(static_cast<std::size_t>(ndof) * static_cast<std::size_t>(ndof), 0.0);
  std::array<std::array<Real, 3>, kMaxNodes> g{};

  for (const GaussPoint& gp : gauss_rule(type)) {
    const Shape s = shape(type, gp.xi, gp.et, gp.ze);
    const Real det = physical_gradients(s, coords, g);
    const Real scale = det * gp.w;

    // B is 6 x ndof. Build DB = D * B (6 x ndof), then Ke += scale * Bᵀ (DB).
    // Column block for node k: rows(xx,yy,zz,xy,xz,yz) x cols(u,v,w).
    auto Bval = [&](int row, int k, int dir) -> Real {
      const Real gx = g[static_cast<std::size_t>(k)][0], gy = g[static_cast<std::size_t>(k)][1], gz = g[static_cast<std::size_t>(k)][2];
      switch (row) {
        case 0: return dir == 0 ? gx : 0.0;                 // xx
        case 1: return dir == 1 ? gy : 0.0;                 // yy
        case 2: return dir == 2 ? gz : 0.0;                 // zz
        case 3: return dir == 0 ? gy : (dir == 1 ? gx : 0.0);  // xy
        case 4: return dir == 0 ? gz : (dir == 2 ? gx : 0.0);  // xz
        default: return dir == 1 ? gz : (dir == 2 ? gy : 0.0); // yz
      }
    };

    for (int a = 0; a < ndof; ++a) {
      const int ka = a / kDofsPerNode, da = a % kDofsPerNode;
      // DB column for dof a: (D * B)[row] = sum_r D[row][r] * B[r][a]
      std::array<Real, kVoigt> DBa{};
      for (int row = 0; row < kVoigt; ++row) {
        Real v = 0.0;
        for (int r = 0; r < kVoigt; ++r) v += D[static_cast<std::size_t>(row)][static_cast<std::size_t>(r)] * Bval(r, ka, da);
        DBa[static_cast<std::size_t>(row)] = v;
      }
      for (int b = a; b < ndof; ++b) {
        const int kb = b / kDofsPerNode, db = b % kDofsPerNode;
        Real v = 0.0;
        for (int row = 0; row < kVoigt; ++row) v += Bval(row, kb, db) * DBa[static_cast<std::size_t>(row)];
        Ke[static_cast<std::size_t>(a) * static_cast<std::size_t>(ndof) + static_cast<std::size_t>(b)] += scale * v;
      }
    }
  }
  // Mirror upper triangle to lower (Ke symmetric).
  for (int a = 0; a < ndof; ++a)
    for (int b = a + 1; b < ndof; ++b)
      Ke[static_cast<std::size_t>(b) * static_cast<std::size_t>(ndof) + static_cast<std::size_t>(a)] =
          Ke[static_cast<std::size_t>(a) * static_cast<std::size_t>(ndof) + static_cast<std::size_t>(b)];
  return Ke;
}

}  // namespace cxpp::fem
