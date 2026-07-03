#include "calculixpp/fem/cavity_radiation.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

#include "calculixpp/fem/element.hpp"

namespace cxpp::fem {
namespace {

constexpr Real kPi = 3.14159265358979323846;

Real dot(const Vec3& a, const Vec3& b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// Direct double-area view-factor integrand contribution between two Gauss points.
// dF = cos_i cos_j / (pi r^2) dA_j, where cos_i = (r_hat . n_i), cos_j = (-r_hat . n_j),
// r pointing from patch i's point to patch j's point. Rays that leave through the back
// of either patch (a negative cosine) do not exchange energy and contribute nothing.
Real pair_kernel(const Vec3& xi, const Vec3& ni, const Vec3& xj, const Vec3& nj) {
  const Vec3 r{xj[0] - xi[0], xj[1] - xi[1], xj[2] - xi[2]};
  const Real r2 = dot(r, r);
  if (r2 <= 1e-30) return 0.0;
  const Real inv = 1.0 / std::sqrt(r2);
  const Vec3 rh{r[0] * inv, r[1] * inv, r[2] * inv};
  const Real ci = dot(rh, ni);         // emitting patch: outward toward j
  const Real cj = -dot(rh, nj);        // receiving patch: inward from i
  if (ci <= 0.0 || cj <= 0.0) return 0.0;
  return ci * cj / (kPi * r2);
}

// The centroid of an element's nodes, used to orient face normals into the cavity.
Vec3 element_centroid(const Mesh& mesh, Index ei) {
  const Element& elem = mesh.elements()[static_cast<std::size_t>(ei)];
  const int n = nodes_per_element(elem.type);
  Vec3 c{0, 0, 0};
  for (int i = 0; i < n; ++i) {
    const Index ni = mesh.node_index(elem.nodes[static_cast<std::size_t>(i)]);
    const Vec3& x = mesh.nodes()[static_cast<std::size_t>(ni)].x;
    for (int d = 0; d < 3; ++d) c[static_cast<std::size_t>(d)] += x[static_cast<std::size_t>(d)];
  }
  for (int d = 0; d < 3; ++d) c[static_cast<std::size_t>(d)] /= static_cast<Real>(n);
  return c;
}

// One radiating patch reduced to geometry, gathered from an element face. The face
// Gauss-point normals from the shared integrator follow the element's node ordering,
// whose sign is topology-dependent; a cavity patch must emit AWAY from its own solid,
// so every normal is flipped to point from the element centroid toward the face (out
// of the solid, into the cavity). Flat faces then share one consistent outward sense.
CavityPatch make_patch(const Mesh& mesh, Index ei, int face, Real emis) {
  CavityPatch p;
  p.emis = emis;
  p.surf = face_surface_integrals(mesh, ei, face);
  // Refined face sampling (sub x sub subcells) for the view-factor double-area
  // quadrature: the 1/r^2 kernel is under-integrated by a single 2x2 rule on opposed
  // patches. sub=6 brings the cube opposite-face view factor within ~1% of Hottel.
  const std::vector<FacePoint> pts = face_gauss_points(mesh, ei, face, 6);
  const Vec3 ecen = element_centroid(mesh, ei);

  // Orientation reference: the face's mean point relative to the element centroid.
  Vec3 fmean{0, 0, 0};
  Real wsum = 0.0;
  for (const FacePoint& fp : pts) {
    wsum += fp.w;
    for (int d = 0; d < 3; ++d) fmean[static_cast<std::size_t>(d)] += fp.w * fp.x[static_cast<std::size_t>(d)];
  }
  if (wsum > 0.0)
    for (int d = 0; d < 3; ++d) fmean[static_cast<std::size_t>(d)] /= wsum;
  const Vec3 outref{fmean[0] - ecen[0], fmean[1] - ecen[1], fmean[2] - ecen[2]};

  Vec3 nsum{0, 0, 0};
  for (const FacePoint& fp0 : pts) {
    FacePoint fp = fp0;
    if (dot(fp.n, outref) < 0.0)  // flip inward normals to point out of the solid
      for (int d = 0; d < 3; ++d) fp.n[static_cast<std::size_t>(d)] = -fp.n[static_cast<std::size_t>(d)];
    p.gp_x.push_back(fp.x);
    p.gp_n.push_back(fp.n);
    p.gp_w.push_back(fp.w);
    p.area += fp.w;
    for (int d = 0; d < 3; ++d) {
      p.centroid[static_cast<std::size_t>(d)] += fp.w * fp.x[static_cast<std::size_t>(d)];
      nsum[static_cast<std::size_t>(d)] += fp.w * fp.n[static_cast<std::size_t>(d)];
    }
  }
  if (p.area > 0.0)
    for (int d = 0; d < 3; ++d) p.centroid[static_cast<std::size_t>(d)] /= p.area;
  // Area-weighted mean normal, renormalized (unit for a flat face).
  const Real nm = std::sqrt(dot(nsum, nsum));
  if (nm > 0.0)
    for (int d = 0; d < 3; ++d) p.normal[static_cast<std::size_t>(d)] = nsum[static_cast<std::size_t>(d)] / nm;
  return p;
}

// Raw view factor F_ij = (1/A_i) ∫_Ai ∫_Aj cos_i cos_j / (pi r^2) dA_j dA_i by direct
// double-area quadrature over the two patches' face Gauss points.
Real view_factor(const CavityPatch& pi, const CavityPatch& pj) {
  if (pi.area <= 0.0) return 0.0;
  Real acc = 0.0;
  for (std::size_t a = 0; a < pi.gp_x.size(); ++a)
    for (std::size_t b = 0; b < pj.gp_x.size(); ++b)
      acc += pi.gp_w[a] * pj.gp_w[b] *
             pair_kernel(pi.gp_x[a], pi.gp_n[a], pj.gp_x[b], pj.gp_n[b]);
  return acc / pi.area;
}

// Solve the dense system A X = B in place for the n x n matrix `A` (row-major) and the
// n x m right-hand-side block `B` (row-major), by Gauss-Jordan elimination with partial
// pivoting. On return `B` holds X (each of the m columns solved). A singular pivot row
// (an isolated patch) is skipped, leaving its solution columns zero. Small dense system
// (n = patch count), so an explicit solver keeps the module dependency-free.
void solve_dense_multi(std::vector<Real>& A, std::vector<Real>& B, int n, int m) {
  const auto aij = [&](int i, int j) -> Real& {
    return A[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) + static_cast<std::size_t>(j)];
  };
  const auto bij = [&](int i, int j) -> Real& {
    return B[static_cast<std::size_t>(i) * static_cast<std::size_t>(m) + static_cast<std::size_t>(j)];
  };
  for (int col = 0; col < n; ++col) {
    int piv = col;
    for (int r = col + 1; r < n; ++r)
      if (std::abs(aij(r, col)) > std::abs(aij(piv, col))) piv = r;
    if (piv != col) {
      for (int c = 0; c < n; ++c) std::swap(aij(col, c), aij(piv, c));
      for (int c = 0; c < m; ++c) std::swap(bij(col, c), bij(piv, c));
    }
    const Real d = aij(col, col);
    if (std::abs(d) < 1e-300) continue;  // isolated patch -> leave its columns zero
    for (int r = 0; r < n; ++r) {
      if (r == col) continue;
      const Real fac = aij(r, col) / d;
      if (fac == 0.0) continue;
      for (int c = col; c < n; ++c) aij(r, c) -= fac * aij(col, c);
      for (int c = 0; c < m; ++c) bij(r, c) -= fac * bij(col, c);
    }
  }
  for (int i = 0; i < n; ++i) {
    const Real d = aij(i, i);
    if (std::abs(d) < 1e-300) continue;
    for (int c = 0; c < m; ++c) bij(i, c) /= d;
  }
}

// Net radiative heat flow leaving one patch and its temperature tangent, given the
// solved radiosities J and the tangent block dJ/de (column 1+k of `aug`). GRAY patch
// (eps < 1): Q_i = A_i eps/(1-eps) (e_i - J_i). BLACK patch (eps == 1): J_i = e_i, so
// Q_i = A_i (e_i - Σ_j F_ij J_j) (the incident irradiation). Fills row i of Q / dQdT.
void patch_heat_flow(const Cavity& cav, int i, int n, const std::vector<Real>& e,
                     const std::vector<Real>& de, const std::vector<Real>& J,
                     const std::vector<Real>& aug, int m, std::vector<Real>& Q,
                     std::vector<Real>& dQdT) {
  const CavityPatch& p = cav.patches[static_cast<std::size_t>(i)];
  const Real Ai = p.area, epsi = p.emis;
  const auto Cik = [&](int r, int k) {
    return aug[static_cast<std::size_t>(r) * static_cast<std::size_t>(m) + 1 + k];
  };
  const auto Fij = [&](int r, int c) {
    return cav.F[static_cast<std::size_t>(r) * static_cast<std::size_t>(n) + static_cast<std::size_t>(c)];
  };
  const auto set = [&](int k, Real v) {
    dQdT[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) + static_cast<std::size_t>(k)] = v;
  };

  if (epsi < 1.0 - 1e-12) {
    const Real g = epsi / (1.0 - epsi);
    Q[static_cast<std::size_t>(i)] = Ai * g * (e[static_cast<std::size_t>(i)] - J[static_cast<std::size_t>(i)]);
    for (int k = 0; k < n; ++k) {
      const Real dJik = Cik(i, k) * de[static_cast<std::size_t>(k)];
      const Real dei = (i == k) ? de[static_cast<std::size_t>(i)] : 0.0;
      set(k, Ai * g * (dei - dJik));
    }
    return;
  }
  Real irr = 0.0;  // Σ_j F_ij J_j
  for (int j = 0; j < n; ++j) irr += Fij(i, j) * J[static_cast<std::size_t>(j)];
  Q[static_cast<std::size_t>(i)] = Ai * (e[static_cast<std::size_t>(i)] - irr);
  for (int k = 0; k < n; ++k) {
    Real dsum = 0.0;  // d(Σ_j F_ij J_j)/dT_k = (Σ_j F_ij C_jk) de_k
    for (int j = 0; j < n; ++j) dsum += Fij(i, j) * Cik(j, k);
    dsum *= de[static_cast<std::size_t>(k)];
    const Real dei = (i == k) ? de[static_cast<std::size_t>(i)] : 0.0;
    set(k, Ai * (dei - dsum));
  }
}

}  // namespace

Cavity build_cavity(const Model& model) {
  const Mesh& mesh = model.mesh;
  Cavity cav;
  for (const Radiate& rd : model.radiates) {
    if (!rd.cavity) continue;
    const Index ei = mesh.element_index(rd.elem_id);
    if (ei < 0) throw std::runtime_error("*RADIATE ...,CR references unknown element");
    cav.patches.push_back(make_patch(mesh, ei, rd.face, rd.emissivity));
  }
  const int n = static_cast<int>(cav.patches.size());
  cav.n = n;
  if (n == 0) return cav;

  cav.F.assign(static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 0.0);
  for (int i = 0; i < n; ++i) {
    Real row = 0.0;
    for (int j = 0; j < n; ++j) {
      if (i == j) continue;  // convex/planar patch: no self-view
      const Real f = view_factor(cav.patches[static_cast<std::size_t>(i)],
                                 cav.patches[static_cast<std::size_t>(j)]);
      cav.F[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) + static_cast<std::size_t>(j)] = f;
      row += f;
    }
    // Enforce the summation rule Σ_j F_ij = 1 (closed enclosure) by rescaling the row.
    // Quadrature and the neglected self/opening view leave the raw row shy of 1; the
    // rescale distributes the deficit proportionally, which is the standard closure.
    if (row > 0.0)
      for (int j = 0; j < n; ++j)
        cav.F[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) + static_cast<std::size_t>(j)] /= row;
  }
  return cav;
}

void cavity_heat_flow(const Cavity& cav, const std::vector<Real>& tabs, Real sigma,
                      std::vector<Real>& Q, std::vector<Real>& dQdT) {
  const int n = cav.n;
  Q.assign(static_cast<std::size_t>(n), 0.0);
  dQdT.assign(static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 0.0);
  if (n == 0) return;

  // Emissive power e_i = sigma T_i^4 and its derivative de_i/dT_i = 4 sigma T_i^3.
  std::vector<Real> e(static_cast<std::size_t>(n)), de(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    const Real t = tabs[static_cast<std::size_t>(i)];
    e[static_cast<std::size_t>(i)] = sigma * t * t * t * t;
    de[static_cast<std::size_t>(i)] = 4.0 * sigma * t * t * t;
  }

  // Gray-body radiosity: (I - (1-eps) F) J = eps e. Solve by dense Gauss elimination
  // (cavities are small — O(patches) faces). Then the net flux leaving patch i is
  //   q_i = eps_i/(1-eps_i) (e_i - J_i)      [W/m^2],  Q_i = A_i q_i.
  // For eps_i == 1 (black) J_i == e_i, and q_i = e_i - Σ_j F_ij J_j directly.
  const auto Fij = [&](int i, int j) {
    return cav.F[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) + static_cast<std::size_t>(j)];
  };

  // Radiosity system M J = eps e with M = I - (1-eps) F. We also want the tangent
  // dJ/de = M^{-1} diag(eps), so we solve the augmented block [ eps*e | eps*I ] in one
  // pass: after the solve, column 0 is J and column 1+k is C_·k = dJ/de_k.
  std::vector<Real> M(static_cast<std::size_t>(n) * static_cast<std::size_t>(n), 0.0);
  const int m = 1 + n;
  std::vector<Real> aug(static_cast<std::size_t>(n) * static_cast<std::size_t>(m), 0.0);
  for (int i = 0; i < n; ++i) {
    const Real epsi = cav.patches[static_cast<std::size_t>(i)].emis;
    const Real refl = 1.0 - epsi;  // 1 - eps_i
    for (int j = 0; j < n; ++j)
      M[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) + static_cast<std::size_t>(j)] =
          (i == j ? 1.0 : 0.0) - refl * Fij(i, j);
    aug[static_cast<std::size_t>(i) * static_cast<std::size_t>(m) + 0] = epsi * e[static_cast<std::size_t>(i)];
    aug[static_cast<std::size_t>(i) * static_cast<std::size_t>(m) + 1 + i] = epsi;  // eps_i on column i
  }
  solve_dense_multi(M, aug, n, m);

  std::vector<Real> J(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i)
    J[static_cast<std::size_t>(i)] = aug[static_cast<std::size_t>(i) * static_cast<std::size_t>(m) + 0];

  for (int i = 0; i < n; ++i)
    patch_heat_flow(cav, i, n, e, de, J, aug, m, Q, dQdT);
}

}  // namespace cxpp::fem
