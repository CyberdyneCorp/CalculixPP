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

// Gauss-Legendre 1D nodes on [-1,1]. 1-point (reduced) and 2-point (full) rules
// are all the hex/wedge families need (shape functions are at most quadratic).
constexpr Real kG2 = 0.577350269189625764;  // 1/sqrt(3)

// Hex tensor-product rules. Weights are the product of 1D Gauss weights (all 1 for
// the 2-point rule). Order follows CalculiX gauss.f: xi fastest, then et, then ze.
//   C3D8  full     : 2x2x2 (8 pts)          -> kGaussHex2
//   C3D8R reduced  : 1 centroid point (w=8) -> kGaussHex1
//   C3D20 full     : 3x3x3 (27 pts)         -> kGaussHex3
//   C3D20R reduced : 2x2x2 (8 pts)          -> kGaussHex2 (reduced vs 3x3x3)
constexpr std::array<GaussPoint, 8> kGaussHex2{{
    {-kG2, -kG2, -kG2, 1.0}, {kG2, -kG2, -kG2, 1.0},
    {-kG2, kG2, -kG2, 1.0},  {kG2, kG2, -kG2, 1.0},
    {-kG2, -kG2, kG2, 1.0},  {kG2, -kG2, kG2, 1.0},
    {-kG2, kG2, kG2, 1.0},   {kG2, kG2, kG2, 1.0},
}};
// C3D8R reduced: single centroid point, weight 2*2*2 = 8 (stabilized by hourglass
// control; see build_hourglass / accumulate_hourglass).
constexpr std::array<GaussPoint, 1> kGaussHex1{{{0.0, 0.0, 0.0, 8.0}}};

// 3-point 1D Gauss-Legendre for the 3x3x3 full C3D20 rule.
constexpr Real kG3 = 0.774596669241483377;  // sqrt(3/5)
constexpr Real kW3a = 5.0 / 9.0;             // outer-node weight
constexpr Real kW3b = 8.0 / 9.0;             // centre-node weight
constexpr std::array<Real, 3> kG3nodes{-kG3, 0.0, kG3};
constexpr std::array<Real, 3> kG3wts{kW3a, kW3b, kW3a};

std::array<GaussPoint, 27> make_hex3() {
  std::array<GaussPoint, 27> r{};
  int idx = 0;
  for (int c = 0; c < 3; ++c)
    for (int b = 0; b < 3; ++b)
      for (int a = 0; a < 3; ++a)
        r[static_cast<std::size_t>(idx++)] = {
            kG3nodes[static_cast<std::size_t>(a)], kG3nodes[static_cast<std::size_t>(b)],
            kG3nodes[static_cast<std::size_t>(c)],
            kG3wts[static_cast<std::size_t>(a)] * kG3wts[static_cast<std::size_t>(b)] *
                kG3wts[static_cast<std::size_t>(c)]};
  return r;
}
const std::array<GaussPoint, 27> kGaussHex3 = make_hex3();

// Wedge rules: triangle rule (in xi,et) tensored with a 1D Gauss rule (in ze).
// C3D6 (linear): 1 triangle point x 2 ze points. Triangle centroid weight 1/2,
// ze weights 1 -> w = 1/2. C3D15 (quadratic): 3 triangle points (weight 1/6 each)
// x 2 ze points -> w = 1/6.
constexpr std::array<GaussPoint, 2> kGaussWedge6{{
    {1.0 / 3.0, 1.0 / 3.0, -kG2, 0.5},
    {1.0 / 3.0, 1.0 / 3.0, kG2, 0.5},
}};
constexpr Real kT = 1.0 / 6.0;
constexpr Real kU = 2.0 / 3.0;
constexpr std::array<GaussPoint, 6> kGaussWedge15{{
    {kT, kT, -kG2, kT}, {kU, kT, -kG2, kT}, {kT, kU, -kG2, kT},
    {kT, kT, kG2, kT},  {kU, kT, kG2, kT},  {kT, kU, kG2, kT},
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
  switch (type) {
    case ElementType::C3D4:
      return kGaussC3D4;
    case ElementType::C3D10:
      return kGaussC3D10;
    case ElementType::C3D8:
      return kGaussHex2;
    case ElementType::C3D20:
      return kGaussHex3;
    case ElementType::C3D8R:
      return kGaussHex1;
    case ElementType::C3D20R:
      return kGaussHex2;
    case ElementType::C3D6:
      return kGaussWedge6;
    case ElementType::C3D15:
      return kGaussWedge15;
  }
  return kGaussC3D4;
}

namespace {

// C3D4 linear tetrahedron (shape4tet.f). Natural coords barycentric, corner 1 at
// the origin.
void shape_c3d4(Real xi, Real et, Real ze, Shape& s) {
  s.n = 4;
  s.N = {1.0 - xi - et - ze, xi, et, ze};
  s.dNdxi[0] = {-1, -1, -1};
  s.dNdxi[1] = {1, 0, 0};
  s.dNdxi[2] = {0, 1, 0};
  s.dNdxi[3] = {0, 0, 1};
}

// C3D10 quadratic tetrahedron (shape10tet.f). a = 1 - xi - et - ze.
void shape_c3d10(Real xi, Real et, Real ze, Shape& s) {
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
}

// C3D8 trilinear hexahedron (shape8h.f). -1<=xi,et,ze<=1. Corner sign pattern
// (sg,sh,sr) with N_k = (1+sg*xi)(1+sh*et)(1+sr*ze)/8.
constexpr int kHexSign[8][3] = {
    {-1, -1, -1}, {1, -1, -1}, {1, 1, -1}, {-1, 1, -1},
    {-1, -1, 1},  {1, -1, 1},  {1, 1, 1},  {-1, 1, 1},
};
void shape_c3d8(Real xi, Real et, Real ze, Shape& s) {
  s.n = 8;
  for (int k = 0; k < 8; ++k) {
    const Real sg = kHexSign[k][0], sh = kHexSign[k][1], sr = kHexSign[k][2];
    const Real g = 1.0 + sg * xi, h = 1.0 + sh * et, r = 1.0 + sr * ze;
    s.N[static_cast<std::size_t>(k)] = g * h * r / 8.0;
    s.dNdxi[static_cast<std::size_t>(k)] = {sg * h * r / 8.0, sh * g * r / 8.0,
                                            sr * g * h / 8.0};
  }
}

// C3D20 serendipity hexahedron (shape20h.f). Nodes 1-8 corners, 9-20 edge
// midsides. Edge midside j connects corner endpoints kMidHex[j] and lies on the
// zero-coordinate of one natural axis (kEdgeAxis[j] is that axis: 0=xi,1=et,2=ze).
constexpr int kMidHex[12][2] = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0},  // 9-12: bottom (ze=-1) edges
    {4, 5}, {5, 6}, {6, 7}, {7, 4},  // 13-16: top (ze=+1) edges
    {0, 4}, {1, 5}, {2, 6}, {3, 7},  // 17-20: vertical (xi,et fixed) edges
};
constexpr int kEdgeAxis[12] = {0, 1, 0, 1, 0, 1, 0, 1, 2, 2, 2, 2};
void shape_c3d20(Real xi, Real et, Real ze, Shape& s) {
  s.n = 20;
  const std::array<Real, 3> p{xi, et, ze};
  // Corner nodes: N = (1+sg*xi)(1+sh*et)(1+sr*ze)(sg*xi+sh*et+sr*ze-2)/8.
  for (int k = 0; k < 8; ++k) {
    const std::array<Real, 3> sgn{static_cast<Real>(kHexSign[k][0]),
                                  static_cast<Real>(kHexSign[k][1]),
                                  static_cast<Real>(kHexSign[k][2])};
    const Real g = 1.0 + sgn[0] * xi, h = 1.0 + sgn[1] * et, r = 1.0 + sgn[2] * ze;
    const Real q = sgn[0] * xi + sgn[1] * et + sgn[2] * ze - 2.0;
    s.N[static_cast<std::size_t>(k)] = g * h * r * q / 8.0;
    const std::array<Real, 3> fac{g, h, r};
    for (int d = 0; d < 3; ++d) {
      // d/dp_d [ g h r q /8 ] = sgn_d*(others)*q/8 + g h r * sgn_d /8.
      Real others = 1.0;
      for (int m = 0; m < 3; ++m)
        if (m != d) others *= fac[static_cast<std::size_t>(m)];
      s.dNdxi[static_cast<std::size_t>(k)][static_cast<std::size_t>(d)] =
          (sgn[static_cast<std::size_t>(d)] * others * q +
           g * h * r * sgn[static_cast<std::size_t>(d)]) /
          8.0;
    }
  }
  // Edge midside nodes: on axis a the shape is (1-p_a^2)/2 times the two
  // in-plane (1+s*p) factors /... Standard serendipity edge function:
  //   N = (1-p_a^2) * prod_{m!=a} (1+s_m p_m) / 4.
  for (int j = 0; j < 12; ++j) {
    const int node = 8 + j;
    const int a = kEdgeAxis[j];
    // In-plane signs come from the endpoints (identical on the two non-axis dirs).
    const int c0 = kMidHex[j][0];
    std::array<Real, 3> sgn{static_cast<Real>(kHexSign[c0][0]),
                            static_cast<Real>(kHexSign[c0][1]),
                            static_cast<Real>(kHexSign[c0][2])};
    const Real pa = p[static_cast<std::size_t>(a)];
    Real plane = 1.0;
    for (int m = 0; m < 3; ++m)
      if (m != a) plane *= 1.0 + sgn[static_cast<std::size_t>(m)] * p[static_cast<std::size_t>(m)];
    s.N[static_cast<std::size_t>(node)] = (1.0 - pa * pa) * plane / 4.0;
    for (int d = 0; d < 3; ++d) {
      if (d == a) {
        s.dNdxi[static_cast<std::size_t>(node)][static_cast<std::size_t>(d)] =
            (-2.0 * pa) * plane / 4.0;
      } else {
        Real other = 1.0;
        for (int m = 0; m < 3; ++m)
          if (m != a && m != d)
            other *= 1.0 + sgn[static_cast<std::size_t>(m)] * p[static_cast<std::size_t>(m)];
        s.dNdxi[static_cast<std::size_t>(node)][static_cast<std::size_t>(d)] =
            (1.0 - pa * pa) * sgn[static_cast<std::size_t>(d)] * other / 4.0;
      }
    }
  }
}

// C3D6 linear wedge (shape6w.f). 0<=xi,et; xi+et<=1; -1<=ze<=1. a = 1-xi-et.
void shape_c3d6(Real xi, Real et, Real ze, Shape& s) {
  s.n = 6;
  const Real a = 1.0 - xi - et;
  const Real m = 0.5 * (1.0 - ze), p = 0.5 * (1.0 + ze);
  s.N = {a * m, xi * m, et * m, a * p, xi * p, et * p};
  s.dNdxi[0] = {-m, -m, -0.5 * a};
  s.dNdxi[1] = {m, 0.0, -0.5 * xi};
  s.dNdxi[2] = {0.0, m, -0.5 * et};
  s.dNdxi[3] = {-p, -p, 0.5 * a};
  s.dNdxi[4] = {p, 0.0, 0.5 * xi};
  s.dNdxi[5] = {0.0, p, 0.5 * et};
}

// C3D15 quadratic wedge (shape15w.f). a = 1-xi-et.
void shape_c3d15(Real xi, Real et, Real ze, Shape& s) {
  s.n = 15;
  const Real a = 1.0 - xi - et;
  const Real omz = 1.0 - ze, opz = 1.0 + ze, omz2 = 1.0 - ze * ze;
  s.N = {-0.5 * a * omz * (2 * xi + 2 * et + ze),
         0.5 * xi * omz * (2 * xi - 2 - ze),
         0.5 * et * omz * (2 * et - 2 - ze),
         -0.5 * a * opz * (2 * xi + 2 * et - ze),
         0.5 * xi * opz * (2 * xi - 2 + ze),
         0.5 * et * opz * (2 * et - 2 + ze),
         2 * xi * a * omz,
         2 * xi * et * omz,
         2 * et * a * omz,
         2 * xi * a * opz,
         2 * xi * et * opz,
         2 * et * a * opz,
         a * omz2,
         xi * omz2,
         et * omz2};
  // xi-derivatives (shape15w.f shp(1,*)).
  const std::array<Real, 15> dxi{
      0.5 * omz * (4 * xi + 4 * et + ze - 2),
      0.5 * omz * (4 * xi - ze - 2),
      0.0,
      0.5 * opz * (4 * xi + 4 * et - ze - 2),
      0.5 * opz * (4 * xi + ze - 2),
      0.0,
      2 * omz * (1 - 2 * xi - et),
      2 * et * omz,
      -2 * et * omz,
      2 * opz * (1 - 2 * xi - et),
      2 * et * opz,
      -2 * et * opz,
      -omz2,
      omz2,
      0.0};
  // et-derivatives (by symmetry xi<->et in the wedge base).
  const std::array<Real, 15> det{
      0.5 * omz * (4 * xi + 4 * et + ze - 2),
      0.0,
      0.5 * omz * (4 * et - ze - 2),
      0.5 * opz * (4 * xi + 4 * et - ze - 2),
      0.0,
      0.5 * opz * (4 * et + ze - 2),
      -2 * xi * omz,
      2 * xi * omz,
      2 * omz * (1 - xi - 2 * et),
      -2 * xi * opz,
      2 * xi * opz,
      2 * opz * (1 - xi - 2 * et),
      -omz2,
      0.0,
      omz2};
  // ze-derivatives (differentiate N above w.r.t. ze).
  std::array<Real, 15> dz{};
  dz[0] = -0.5 * a * (-1.0 * (2 * xi + 2 * et + ze) + omz * 1.0);
  dz[1] = 0.5 * xi * (-1.0 * (2 * xi - 2 - ze) + omz * (-1.0));
  dz[2] = 0.5 * et * (-1.0 * (2 * et - 2 - ze) + omz * (-1.0));
  dz[3] = -0.5 * a * (1.0 * (2 * xi + 2 * et - ze) + opz * (-1.0));
  dz[4] = 0.5 * xi * (1.0 * (2 * xi - 2 + ze) + opz * (1.0));
  dz[5] = 0.5 * et * (1.0 * (2 * et - 2 + ze) + opz * (1.0));
  dz[6] = -2 * xi * a;
  dz[7] = -2 * xi * et;
  dz[8] = -2 * et * a;
  dz[9] = 2 * xi * a;
  dz[10] = 2 * xi * et;
  dz[11] = 2 * et * a;
  dz[12] = a * (-2 * ze);
  dz[13] = xi * (-2 * ze);
  dz[14] = et * (-2 * ze);
  for (int k = 0; k < 15; ++k)
    s.dNdxi[static_cast<std::size_t>(k)] = {dxi[static_cast<std::size_t>(k)],
                                            det[static_cast<std::size_t>(k)],
                                            dz[static_cast<std::size_t>(k)]};
}

}  // namespace

Shape shape(ElementType type, Real xi, Real et, Real ze) {
  Shape s;
  switch (type) {
    case ElementType::C3D4:
      shape_c3d4(xi, et, ze, s);
      break;
    case ElementType::C3D10:
      shape_c3d10(xi, et, ze, s);
      break;
    case ElementType::C3D8:
    case ElementType::C3D8R:
      shape_c3d8(xi, et, ze, s);
      break;
    case ElementType::C3D20:
    case ElementType::C3D20R:
      shape_c3d20(xi, et, ze, s);
      break;
    case ElementType::C3D6:
      shape_c3d6(xi, et, ze, s);
      break;
    case ElementType::C3D15:
      shape_c3d15(xi, et, ze, s);
      break;
  }
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

Real b_entry(int row, const std::array<Real, 3>& gk, int dir) {
  const Real gx = gk[0], gy = gk[1], gz = gk[2];
  switch (row) {
    case 0: return dir == 0 ? gx : 0.0;                     // xx
    case 1: return dir == 1 ? gy : 0.0;                     // yy
    case 2: return dir == 2 ? gz : 0.0;                     // zz
    case 3: return dir == 0 ? gy : (dir == 1 ? gx : 0.0);   // xy
    case 4: return dir == 0 ? gz : (dir == 2 ? gx : 0.0);   // xz
    default: return dir == 1 ? gz : (dir == 2 ? gy : 0.0);  // yz
  }
}

Voigt6 strain_from_gradients(const std::array<std::array<Real, 3>, kMaxNodes>& g,
                             int n, std::span<const Vec3> ue) {
  // du[i][j] = sum_k g_k[j] * u_k[i]; engineering shear.
  std::array<std::array<Real, 3>, 3> du{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      Real v = 0.0;
      for (int k = 0; k < n; ++k)
        v += g[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)] *
             ue[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)];
      du[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = v;
    }
  return {du[0][0], du[1][1], du[2][2], du[0][1] + du[1][0],
          du[0][2] + du[2][0], du[1][2] + du[2][1]};
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

namespace {

// Accumulate one Gauss point's B^T D B (upper triangle) into Ke, exactly as the
// original element_stiffness did: DB = D*B per column, then Ke[a][b] += scale*(Bᵀ DB)
// for b >= a. Kept bit-identical so element_stiffness results are unchanged.
void accumulate_ktangent(const std::array<std::array<Real, 3>, kMaxNodes>& g,
                         int ndof, const D6& D, Real scale, std::vector<Real>& Ke) {
  for (int a = 0; a < ndof; ++a) {
    const int ka = a / kDofsPerNode, da = a % kDofsPerNode;
    std::array<Real, kVoigt> DBa{};
    for (int row = 0; row < kVoigt; ++row) {
      Real v = 0.0;
      for (int r = 0; r < kVoigt; ++r)
        v += D[static_cast<std::size_t>(row)][static_cast<std::size_t>(r)] *
             b_entry(r, g[static_cast<std::size_t>(ka)], da);
      DBa[static_cast<std::size_t>(row)] = v;
    }
    for (int b = a; b < ndof; ++b) {
      const int kb = b / kDofsPerNode, db = b % kDofsPerNode;
      Real v = 0.0;
      for (int row = 0; row < kVoigt; ++row)
        v += b_entry(row, g[static_cast<std::size_t>(kb)], db) *
             DBa[static_cast<std::size_t>(row)];
      Ke[static_cast<std::size_t>(a) * static_cast<std::size_t>(ndof) +
         static_cast<std::size_t>(b)] += scale * v;
    }
  }
}

// Accumulate one Gauss point's B^T sigma into fe (size ndof).
void accumulate_fint(const std::array<std::array<Real, 3>, kMaxNodes>& g, int ndof,
                     const Voigt6& sigma, Real scale, std::vector<Real>& fe) {
  for (int a = 0; a < ndof; ++a) {
    const int ka = a / kDofsPerNode, da = a % kDofsPerNode;
    Real v = 0.0;
    for (int row = 0; row < kVoigt; ++row)
      v += b_entry(row, g[static_cast<std::size_t>(ka)], da) *
           sigma[static_cast<std::size_t>(row)];
    fe[static_cast<std::size_t>(a)] += scale * v;
  }
}

// Mirror the upper triangle of an ndof x ndof row-major matrix to the lower.
void mirror_upper(int ndof, std::vector<Real>& Ke) {
  for (int a = 0; a < ndof; ++a)
    for (int b = a + 1; b < ndof; ++b)
      Ke[static_cast<std::size_t>(b) * static_cast<std::size_t>(ndof) +
         static_cast<std::size_t>(a)] =
          Ke[static_cast<std::size_t>(a) * static_cast<std::size_t>(ndof) +
             static_cast<std::size_t>(b)];
}

// ---- Hourglass control for reduced-integration C3D8R --------------------------
//
// A single-point-integrated C3D8 stiffness is rank-deficient: beyond the 6
// rigid-body modes it admits 12 spurious zero-energy "hourglass" modes (4 per
// spatial direction) that produce no strain at the centroid. We stabilize them
// with Flanagan-Belytschko hourglass control (Flanagan & Belytschko, IJNME 1981).
//
// The 4 hourglass base vectors Γ_α (length-8, per corner) span the space
// orthogonal to the constant and linear fields at the centroid:
constexpr Real kHgGamma[4][8] = {
    {1, 1, -1, -1, -1, -1, 1, 1},   // h1
    {1, -1, -1, 1, -1, 1, 1, -1},   // h2
    {1, -1, 1, -1, 1, -1, 1, -1},   // h3
    {-1, 1, -1, 1, 1, -1, 1, -1},   // h4
};

// Build the corrected hourglass shape vectors γ_α (Flanagan-Belytschko), the
// centroid B-matrix b_j (mean gradients), the element volume, and a scalar
// stabilization stiffness κ. Only meaningful for an 8-node hex. Returns false if
// the element is degenerate (non-positive centroid Jacobian, already thrown).
struct HourglassData {
  Real gamma[4][8];  // corrected hourglass vectors
  Real kappa;        // stabilization modulus (per hourglass mode, per direction)
};

HourglassData build_hourglass(std::span<const Vec3> coords, Real mu) {
  // Centroid gradients b_j = dN_j/dX at (0,0,0) and element volume.
  const Shape s = shape(ElementType::C3D8R, 0.0, 0.0, 0.0);
  std::array<std::array<Real, 3>, kMaxNodes> b{};
  const Real det = physical_gradients(s, coords, b);
  const Real vol = 8.0 * det;  // reduced rule weight is 8

  HourglassData hg;
  for (int al = 0; al < 4; ++al) {
    // γ_α = Γ_α - (Γ_α · x_i) b_i   (removes the part of Γ representable by the
    // linear field, so γ is orthogonal to rigid-body + uniform-strain modes).
    Real gx[3] = {0, 0, 0};
    for (int j = 0; j < 8; ++j)
      for (int d = 0; d < 3; ++d)
        gx[d] += kHgGamma[al][j] * coords[static_cast<std::size_t>(j)][static_cast<std::size_t>(d)];
    for (int j = 0; j < 8; ++j) {
      Real corr = 0.0;
      for (int d = 0; d < 3; ++d) corr += gx[d] * b[static_cast<std::size_t>(j)][static_cast<std::size_t>(d)];
      hg.gamma[al][j] = kHgGamma[al][j] - corr;
    }
  }

  // Stabilization modulus. A physically-scaled choice (Belytschko): κ scales with
  // the shear modulus, the element volume, and the squared mean gradient. We use
  // κ = c_hg * μ * V * Σ_j |b_j|^2 with a small dimensionless c_hg so the added
  // stiffness suppresses hourglassing without over-stiffening the physical modes.
  Real b2 = 0.0;
  for (int j = 0; j < 8; ++j)
    for (int d = 0; d < 3; ++d)
      b2 += b[static_cast<std::size_t>(j)][static_cast<std::size_t>(d)] *
            b[static_cast<std::size_t>(j)][static_cast<std::size_t>(d)];
  constexpr Real c_hg = 0.05;  // documented hourglass coefficient (see notes)
  hg.kappa = c_hg * mu * vol * b2;
  return hg;
}

// Add the hourglass stabilization stiffness κ Σ_α (γ_α ⊗ γ_α) ⊗ I_3 into the upper
// triangle of the 24x24 C3D8R element stiffness.
void accumulate_hourglass(const HourglassData& hg, std::vector<Real>& Ke) {
  constexpr int ndof = 24;
  for (int al = 0; al < 4; ++al)
    for (int i = 0; i < 8; ++i)
      for (int j = 0; j < 8; ++j) {
        const Real kij = hg.kappa * hg.gamma[al][i] * hg.gamma[al][j];
        for (int d = 0; d < 3; ++d) {
          const int a = i * 3 + d, bcol = j * 3 + d;
          if (bcol >= a)
            Ke[static_cast<std::size_t>(a) * ndof + static_cast<std::size_t>(bcol)] += kij;
        }
      }
}

// Hourglass resisting force f_hg = K_hg u for the material-point (internal-force)
// path. K_hg u = κ Σ_α γ_α (γ_α · u_d) per direction d.
void accumulate_hourglass_force(const HourglassData& hg, std::span<const Vec3> ue,
                                std::vector<Real>& fe) {
  for (int al = 0; al < 4; ++al) {
    Real gu[3] = {0, 0, 0};
    for (int j = 0; j < 8; ++j)
      for (int d = 0; d < 3; ++d)
        gu[d] += hg.gamma[al][j] * ue[static_cast<std::size_t>(j)][static_cast<std::size_t>(d)];
    for (int i = 0; i < 8; ++i)
      for (int d = 0; d < 3; ++d)
        fe[static_cast<std::size_t>(i * 3 + d)] += hg.kappa * hg.gamma[al][i] * gu[d];
  }
}

}  // namespace

std::vector<Real> element_stiffness(ElementType type, std::span<const Vec3> coords,
                                    const ElasticIso& mat) {
  const int n = nodes_per_element(type);
  const int ndof = n * kDofsPerNode;
  const D6 D = elastic_iso_D(mat);
  std::vector<Real> Ke(
      static_cast<std::size_t>(ndof) * static_cast<std::size_t>(ndof), 0.0);
  std::array<std::array<Real, 3>, kMaxNodes> g{};

  for (const GaussPoint& gp : gauss_rule(type)) {
    const Shape s = shape(type, gp.xi, gp.et, gp.ze);
    const Real det = physical_gradients(s, coords, g);
    accumulate_ktangent(g, ndof, D, det * gp.w, Ke);
  }
  if (type == ElementType::C3D8R) {
    const Real mu = mat.E / (2.0 * (1.0 + mat.nu));
    accumulate_hourglass(build_hourglass(coords, mu), Ke);
  }
  mirror_upper(ndof, Ke);
  return Ke;
}

std::vector<Real> element_thermal_load(ElementType type,
                                       std::span<const Vec3> coords, const D6& D,
                                       Real alpha, std::span<const Real> te) {
  const int n = nodes_per_element(type);
  const int ndof = n * kDofsPerNode;
  std::vector<Real> fe(static_cast<std::size_t>(ndof), 0.0);
  std::array<std::array<Real, 3>, kMaxNodes> g{};

  for (const GaussPoint& gp : gauss_rule(type)) {
    const Shape s = shape(type, gp.xi, gp.et, gp.ze);
    const Real det = physical_gradients(s, coords, g);
    // Temperature change interpolated to this Gauss point: dT = Σ_k N_k te_k.
    Real dT = 0.0;
    for (int k = 0; k < n; ++k)
      dT += s.N[static_cast<std::size_t>(k)] * te[static_cast<std::size_t>(k)];
    // Thermal strain (normal components only) and its stress D eps_th.
    const Real eth = alpha * dT;
    const Voigt6 eps_th{eth, eth, eth, 0.0, 0.0, 0.0};
    Voigt6 sig_th{};
    for (int i = 0; i < kVoigt; ++i) {
      Real v = 0.0;
      for (int j = 0; j < kVoigt; ++j)
        v += D[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] *
             eps_th[static_cast<std::size_t>(j)];
      sig_th[static_cast<std::size_t>(i)] = v;
    }
    accumulate_fint(g, ndof, sig_th, det * gp.w, fe);
  }
  return fe;
}

ElementResponse element_tangent_force(ElementType type, std::span<const Vec3> coords,
                                      std::span<const Vec3> ue,
                                      const MaterialModel& material,
                                      std::span<MaterialState> state) {
  const int n = nodes_per_element(type);
  const int ndof = n * kDofsPerNode;
  ElementResponse out;
  out.Ke.assign(static_cast<std::size_t>(ndof) * static_cast<std::size_t>(ndof), 0.0);
  out.fe.assign(static_cast<std::size_t>(ndof), 0.0);
  std::array<std::array<Real, 3>, kMaxNodes> g{};

  const auto rule = gauss_rule(type);
  Real mu_hg = 0.0;  // shear modulus from the tangent, for C3D8R hourglass control
  for (std::size_t q = 0; q < rule.size(); ++q) {
    const GaussPoint& gp = rule[q];
    const Shape s = shape(type, gp.xi, gp.et, gp.ze);
    const Real det = physical_gradients(s, coords, g);
    const Real scale = det * gp.w;

    const Voigt6 strain = strain_from_gradients(g, n, ue);
    const MaterialResponse resp = material.evaluate(strain, state[q]);
    accumulate_ktangent(g, ndof, resp.tangent, scale, out.Ke);
    accumulate_fint(g, ndof, resp.stress, scale, out.fe);
    mu_hg = resp.tangent[3][3];  // D[xy][xy] = shear modulus (isotropic)
  }
  // C3D8R hourglass stabilization: add K_hg to the tangent and the resisting force
  // K_hg u to the internal force, so the reduced element is rank-sufficient and the
  // internal force stays consistent (fe = (K_gp + K_hg) ue for linear elasticity).
  if (type == ElementType::C3D8R) {
    const HourglassData hg = build_hourglass(coords, mu_hg);
    accumulate_hourglass(hg, out.Ke);
    accumulate_hourglass_force(hg, ue, out.fe);
  }
  mirror_upper(ndof, out.Ke);
  return out;
}

std::vector<Real> element_conduction(ElementType type, std::span<const Vec3> coords,
                                     Real k) {
  const int n = nodes_per_element(type);
  std::vector<Real> Kt(static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 0.0);
  std::array<std::array<Real, 3>, kMaxNodes> g{};
  for (const GaussPoint& gp : gauss_rule(type)) {
    const Shape s = shape(type, gp.xi, gp.et, gp.ze);
    const Real det = physical_gradients(s, coords, g);
    const Real scale = k * det * gp.w;
    // Kt[a][b] += (g_a · g_b) * k * detJ * w. Symmetric; fill the full matrix.
    for (int a = 0; a < n; ++a)
      for (int b = 0; b < n; ++b) {
        const Real dot = g[static_cast<std::size_t>(a)][0] * g[static_cast<std::size_t>(b)][0] +
                         g[static_cast<std::size_t>(a)][1] * g[static_cast<std::size_t>(b)][1] +
                         g[static_cast<std::size_t>(a)][2] * g[static_cast<std::size_t>(b)][2];
        Kt[static_cast<std::size_t>(a) * static_cast<std::size_t>(n) +
           static_cast<std::size_t>(b)] += scale * dot;
      }
  }
  return Kt;
}

std::vector<Real> element_capacitance(ElementType type, std::span<const Vec3> coords,
                                      Real rho_c) {
  const int n = nodes_per_element(type);
  std::vector<Real> Ce(static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 0.0);
  if (rho_c == 0.0) return Ce;
  std::array<std::array<Real, 3>, kMaxNodes> g{};
  for (const GaussPoint& gp : gauss_rule(type)) {
    const Shape s = shape(type, gp.xi, gp.et, gp.ze);
    const Real det = physical_gradients(s, coords, g);  // detJ only
    const Real scale = rho_c * det * gp.w;
    for (int a = 0; a < n; ++a)
      for (int b = 0; b < n; ++b)
        Ce[static_cast<std::size_t>(a) * static_cast<std::size_t>(n) +
           static_cast<std::size_t>(b)] +=
            scale * s.N[static_cast<std::size_t>(a)] * s.N[static_cast<std::size_t>(b)];
  }
  return Ce;
}

ElasticIsoMaterial::ElasticIsoMaterial(const ElasticIso& props)
    : D_(elastic_iso_D(props)) {}

MaterialResponse ElasticIsoMaterial::evaluate(const Voigt6& strain,
                                              MaterialState& /*state*/) const {
  MaterialResponse r;
  r.tangent = D_;
  for (int i = 0; i < kVoigt; ++i) {
    Real v = 0.0;
    for (int j = 0; j < kVoigt; ++j)
      v += D_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] *
           strain[static_cast<std::size_t>(j)];
    r.stress[static_cast<std::size_t>(i)] = v;
  }
  return r;
}

}  // namespace cxpp::fem
