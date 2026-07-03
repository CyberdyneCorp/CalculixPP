#include "calculixpp/fem/loads.hpp"

#include <array>
#include <cmath>
#include <span>
#include <stdexcept>

#include "calculixpp/fem/element.hpp"

namespace cxpp::fem {
namespace {

// Tet face node ordering (0-based), from CalculiX ifacet (1-based:
// F1 1,3,2,7,6,5 / F2 1,2,4,5,9,8 / F3 2,3,4,6,10,9 / F4 1,4,3,8,10,7).
// Corners are ordered so that (b-a)x(c-a) points outward.
constexpr int kFace[4][6] = {
    {0, 2, 1, 6, 5, 4},
    {0, 1, 3, 4, 8, 7},
    {1, 2, 3, 5, 9, 8},
    {0, 3, 2, 7, 9, 6},
};

// Hex face node ordering (0-based), CalculiX ifaceq (1-based, corners CCW so the
// (t_xi x t_eta) normal points OUTWARD, then edge midsides for C3D20):
//   F1 1,2,3,4 (9,10,11,12)   F2 5,8,7,6 (16,15,14,13)
//   F3 1,5,6,2 (17,13,18,9)   F4 2,6,7,3 (18,14,19,10)
//   F5 3,7,8,4 (19,15,20,11)  F6 4,8,5,1 (20,16,17,12)
constexpr int kFaceHex[6][8] = {
    {0, 1, 2, 3, 8, 9, 10, 11},
    {4, 7, 6, 5, 15, 14, 13, 12},
    {0, 4, 5, 1, 16, 12, 17, 8},
    {1, 5, 6, 2, 17, 13, 18, 9},
    {2, 6, 7, 3, 18, 14, 19, 10},
    {3, 7, 4, 0, 19, 15, 16, 11},
};

// Wedge face node ordering (0-based), CalculiX ifacew (1-based; two triangular end
// faces F1,F2 and three quad side faces F3,F4,F5). Midside nodes appended for C3D15:
//   F1 1,3,2 (9,8,7)          F2 4,5,6 (13,14,15)
//   F3 1,2,5,4 (7,11,13,10)   F4 2,3,6,5 (8,12,14,11)  F5 3,1,4,6 (9,10,15,12)
constexpr int kFaceWedge[5][8] = {
    {0, 2, 1, 8, 7, 6, -1, -1},
    {3, 4, 5, 12, 13, 14, -1, -1},
    {0, 1, 4, 3, 6, 10, 12, 9},
    {1, 2, 5, 4, 7, 11, 13, 10},
    {2, 0, 3, 5, 8, 9, 14, 11},
};

// A face quadrature point in the face's natural (xi,eta) coordinates.
struct FaceGP {
  Real xi, eta, w;
};

// 3-point order-2 triangle Gauss rule (reference triangle area 1/2). Used for the
// triangular faces of tets (T3/T6) and wedge end faces.
constexpr std::array<FaceGP, 3> kTri{{
    {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0},
    {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
    {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
}};

// 2x2 Gauss rule on the [-1,1]^2 reference quad for hex/wedge quadrilateral faces.
constexpr Real kQg = 0.577350269189625764;  // 1/sqrt(3)
constexpr std::array<FaceGP, 4> kQuad{{
    {-kQg, -kQg, 1.0}, {kQg, -kQg, 1.0}, {kQg, kQg, 1.0}, {-kQg, kQg, 1.0},
}};

// The face geometry a load card sees: how many face nodes, whether the face is a
// triangle or a quad, and which quadrature to use. Triangular faces use barycentric
// (xi,eta) on the unit triangle; quad faces use (xi,eta) on [-1,1]^2.
struct FaceKind {
  bool quad;  // true -> Q4/Q8 quad face; false -> T3/T6 triangle face
  int nf;     // node count (3/6 triangle, 4/8 quad)
};

// Face shape functions and derivatives for a triangle (T3 nf=3, T6 nf=6) or quad
// (Q4 nf=4, Q8 nf=8) at natural coords (a,b). Fills up to 8 entries.
void face_shape(const FaceKind& fk, Real a, Real b, std::array<Real, 8>& N,
                std::array<Real, 8>& dxi, std::array<Real, 8>& deta) {
  N = {}; dxi = {}; deta = {};
  if (!fk.quad) {
    const Real l1 = 1.0 - a - b, l2 = a, l3 = b;
    if (fk.nf == 3) {
      N = {l1, l2, l3, 0, 0, 0, 0, 0};
      dxi = {-1, 1, 0, 0, 0, 0, 0, 0};
      deta = {-1, 0, 1, 0, 0, 0, 0, 0};
      return;
    }
    N = {l1 * (2 * l1 - 1), l2 * (2 * l2 - 1), l3 * (2 * l3 - 1),
         4 * l1 * l2,       4 * l2 * l3,       4 * l3 * l1, 0, 0};
    dxi = {-(4 * l1 - 1), 4 * l2 - 1, 0, 4 * (l1 - l2), 4 * l3, -4 * l3, 0, 0};
    deta = {-(4 * l1 - 1), 0, 4 * l3 - 1, -4 * l2, 4 * l2, 4 * (l1 - l3), 0, 0};
    return;
  }
  // Quad faces on [-1,1]^2. Q4 bilinear corners; Q8 serendipity (corners + 4 edges).
  constexpr int sgn[4][2] = {{-1, -1}, {1, -1}, {1, 1}, {-1, 1}};
  if (fk.nf == 4) {
    for (int k = 0; k < 4; ++k) {
      const Real g = 1.0 + sgn[k][0] * a, h = 1.0 + sgn[k][1] * b;
      N[static_cast<std::size_t>(k)] = 0.25 * g * h;
      dxi[static_cast<std::size_t>(k)] = 0.25 * sgn[k][0] * h;
      deta[static_cast<std::size_t>(k)] = 0.25 * sgn[k][1] * g;
    }
    return;
  }
  // Q8 serendipity: corners 0-3, edge midsides 4-7 (edges 0-1,1-2,2-3,3-0).
  for (int k = 0; k < 4; ++k) {
    const Real sa = sgn[k][0], sb = sgn[k][1];
    const Real g = 1.0 + sa * a, h = 1.0 + sb * b;
    N[static_cast<std::size_t>(k)] = 0.25 * g * h * (sa * a + sb * b - 1.0);
    dxi[static_cast<std::size_t>(k)] = 0.25 * (sa * h * (sa * a + sb * b - 1.0) + g * h * sa);
    deta[static_cast<std::size_t>(k)] = 0.25 * (sb * g * (sa * a + sb * b - 1.0) + g * h * sb);
  }
  // Edge midsides: 4 on eta-const edges (xi free) alternate with xi-const edges.
  constexpr int edgeAxis[4] = {0, 1, 0, 1};  // node 4:xi-edge(b=-1),5:eta-edge(a=1),6:xi(b=1),7:eta(a=-1)
  constexpr Real edgeFix[4] = {-1, 1, 1, -1};
  for (int e = 0; e < 4; ++e) {
    const int node = 4 + e;
    if (edgeAxis[e] == 0) {  // midside on b = edgeFix (xi varies): N = (1-a^2)(1+sb*b)/2
      const Real sb = edgeFix[e];
      N[static_cast<std::size_t>(node)] = 0.5 * (1.0 - a * a) * (1.0 + sb * b);
      dxi[static_cast<std::size_t>(node)] = 0.5 * (-2.0 * a) * (1.0 + sb * b);
      deta[static_cast<std::size_t>(node)] = 0.5 * (1.0 - a * a) * sb;
    } else {  // midside on a = edgeFix (eta varies): N = (1+sa*a)(1-b^2)/2
      const Real sa = edgeFix[e];
      N[static_cast<std::size_t>(node)] = 0.5 * (1.0 + sa * a) * (1.0 - b * b);
      dxi[static_cast<std::size_t>(node)] = 0.5 * sa * (1.0 - b * b);
      deta[static_cast<std::size_t>(node)] = 0.5 * (1.0 + sa * a) * (-2.0 * b);
    }
  }
}

// Resolve the face topology (node list + kind) for one element face. Supports tet
// triangular faces, hex quad faces, and wedge (2 triangle + 3 quad) faces.
FaceKind face_kind_of(ElementType type, int face, std::vector<int>& fn);

// Face geometry gathered for a surface-integral load: local->global node indices
// and the face node coordinates. Shared by pressure and flux integration.
struct FaceGeom {
  std::vector<int> fn;       // element-local face node indices
  std::vector<Index> gnode;  // global (mesh) node indices
  std::vector<Vec3> x;       // face node coordinates
  FaceKind kind;             // triangle/quad + node count
};

FaceGeom gather_face(const Mesh& mesh, const Element& elem, int face) {
  FaceGeom fg;
  fg.kind = face_kind_of(elem.type, face, fg.fn);
  const int nf = fg.kind.nf;
  fg.gnode.resize(static_cast<std::size_t>(nf));
  fg.x.resize(static_cast<std::size_t>(nf));
  for (int i = 0; i < nf; ++i) {
    const Index ni = mesh.node_index(
        elem.nodes[static_cast<std::size_t>(fg.fn[static_cast<std::size_t>(i)])]);
    fg.gnode[static_cast<std::size_t>(i)] = ni;
    fg.x[static_cast<std::size_t>(i)] = mesh.nodes()[static_cast<std::size_t>(ni)].x;
  }
  return fg;
}

// The face quadrature for a face kind (triangle -> kTri, quad -> kQuad).
std::span<const FaceGP> face_gauss(const FaceKind& fk) {
  if (fk.quad) return std::span<const FaceGP>(kQuad.data(), kQuad.size());
  return std::span<const FaceGP>(kTri.data(), kTri.size());
}

// At one face quadrature point, evaluate the shape N and the covariant tangents
// t_xi, t_eta so callers can form the outward normal (pressure) or |dA| (flux).
void face_point(const FaceGeom& fg, const FaceGP& gp, std::array<Real, 8>& N,
                Vec3& txi, Vec3& teta) {
  std::array<Real, 8> dxi{}, deta{};
  face_shape(fg.kind, gp.xi, gp.eta, N, dxi, deta);
  txi = {0, 0, 0};
  teta = {0, 0, 0};
  for (int i = 0; i < fg.kind.nf; ++i)
    for (int d = 0; d < 3; ++d) {
      txi[static_cast<std::size_t>(d)] += dxi[static_cast<std::size_t>(i)] *
                                          fg.x[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)];
      teta[static_cast<std::size_t>(d)] += deta[static_cast<std::size_t>(i)] *
                                           fg.x[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)];
    }
}

// Accumulate consistent nodal forces for a pressure `p` on element face (*DLOAD
// P<face> or *DSLOAD), F_i = -p * integral( N_i * (t_xi x t_eta) ), into f.
void accumulate_pressure_face(const Mesh& mesh, Index ei, int face, Real p,
                              std::vector<Real>& f) {
  const Element& elem = mesh.elements()[static_cast<std::size_t>(ei)];
  const FaceGeom fg = gather_face(mesh, elem, face);
  const int nf = fg.kind.nf;

  std::array<Real, 8> N{};
  Vec3 txi{}, teta{};
  std::vector<Vec3> force(static_cast<std::size_t>(nf), Vec3{0, 0, 0});
  for (const FaceGP& gp : face_gauss(fg.kind)) {
    face_point(fg, gp, N, txi, teta);
    const Vec3 nrm{txi[1] * teta[2] - txi[2] * teta[1],
                   txi[2] * teta[0] - txi[0] * teta[2],
                   txi[0] * teta[1] - txi[1] * teta[0]};
    for (int i = 0; i < nf; ++i)
      for (int d = 0; d < 3; ++d)
        force[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)] +=
            -p * N[static_cast<std::size_t>(i)] * nrm[static_cast<std::size_t>(d)] * gp.w;
  }

  for (int i = 0; i < nf; ++i)
    for (int d = 0; d < 3; ++d)
      f[static_cast<std::size_t>(fg.gnode[static_cast<std::size_t>(i)]) * kDofsPerNode + static_cast<std::size_t>(d)] +=
          force[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)];
}

// Accumulate consistent nodal heat flux for a distributed surface flux `flux`
// (*DFLUX S<face>) into the scalar thermal rhs q (size num_nodes). Positive flux
// flows into the element: q_i += flux * integral( N_i dA ) = flux * Σ_gp N_i |dA| w,
// with |dA| = |t_xi x t_eta| the surface Jacobian. (spec: heat-transfer-analysis.)
void accumulate_flux_face(const Mesh& mesh, Index ei, int face, Real flux,
                          std::vector<Real>& q) {
  const Element& elem = mesh.elements()[static_cast<std::size_t>(ei)];
  const FaceGeom fg = gather_face(mesh, elem, face);
  const int nf = fg.kind.nf;

  std::array<Real, 8> N{};
  Vec3 txi{}, teta{};
  for (const FaceGP& gp : face_gauss(fg.kind)) {
    face_point(fg, gp, N, txi, teta);
    const Vec3 nrm{txi[1] * teta[2] - txi[2] * teta[1],
                   txi[2] * teta[0] - txi[0] * teta[2],
                   txi[0] * teta[1] - txi[1] * teta[0]};
    const Real dA = std::sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]);
    for (int i = 0; i < nf; ++i)
      q[static_cast<std::size_t>(fg.gnode[static_cast<std::size_t>(i)])] +=
          flux * N[static_cast<std::size_t>(i)] * dA * gp.w;
  }
}

Vec3 normalize(Vec3 v) {
  const Real m = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  if (m <= 0.0) return v;
  return {v[0] / m, v[1] / m, v[2] / m};
}

// Body-force density per unit volume at a physical point `x` for one body load:
// GRAV  -> rho * g * dir_hat  (constant over the element).
// CENTRIF -> rho * omega^2 * r_perp, where r_perp is the vector from the axis to
// the point (outward radial direction), integrated with the element shape rule.
Vec3 body_force_density(const BodyLoad& bl, Real rho, Real factor, const Vec3& x) {
  if (bl.kind == BodyLoad::Kind::Gravity) {
    const Vec3 d = normalize(bl.dir);
    const Real s = rho * factor * bl.magnitude;
    return {s * d[0], s * d[1], s * d[2]};
  }
  // Centrifugal: r_perp = (x - point) - ((x - point) . axis_hat) axis_hat.
  const Vec3 axis = normalize(bl.dir);
  const Vec3 rel{x[0] - bl.point[0], x[1] - bl.point[1], x[2] - bl.point[2]};
  const Real along = rel[0] * axis[0] + rel[1] * axis[1] + rel[2] * axis[2];
  const Vec3 rperp{rel[0] - along * axis[0], rel[1] - along * axis[1],
                   rel[2] - along * axis[2]};
  const Real s = rho * factor * bl.magnitude;  // magnitude == omega^2
  return {s * rperp[0], s * rperp[1], s * rperp[2]};
}

// Add one body load's consistent nodal forces over its element set to f.
void accumulate_body_load(const Model& model, const BodyLoad& bl, Real factor,
                          std::vector<Real>& f) {
  const Mesh& mesh = model.mesh;
  const std::vector<Real> rho = model.element_density();
  const std::vector<Index>* elset =
      bl.elset.empty() ? nullptr : mesh.elset(bl.elset);
  if (!bl.elset.empty() && elset == nullptr)
    throw std::runtime_error("body load references unknown elset '" + bl.elset + "'");

  const std::vector<bool> active = model.element_active_mask();
  std::array<std::array<Real, 3>, kMaxNodes> grad{};
  const auto do_element = [&](std::size_t e) {
    if (!active[e]) return;      // *MODEL CHANGE, REMOVE: no body load
    if (rho[e] <= 0.0) return;  // no *DENSITY -> no body force
    const Element& elem = mesh.elements()[e];
    const int n = nodes_per_element(elem.type);
    std::vector<Vec3> coords(static_cast<std::size_t>(n));
    std::vector<Index> nidx(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
      const Index ni = mesh.node_index(elem.nodes[static_cast<std::size_t>(i)]);
      nidx[static_cast<std::size_t>(i)] = ni;
      coords[static_cast<std::size_t>(i)] = mesh.nodes()[static_cast<std::size_t>(ni)].x;
    }
    for (const GaussPoint& gp : gauss_rule(elem.type)) {
      const Shape s = shape(elem.type, gp.xi, gp.et, gp.ze);
      const Real det = physical_gradients(s, coords, grad);
      const Real dv = det * gp.w;
      // Physical point at the Gauss point (for centrifugal radius).
      Vec3 xg{0, 0, 0};
      for (int i = 0; i < n; ++i)
        for (int d = 0; d < 3; ++d)
          xg[static_cast<std::size_t>(d)] += s.N[static_cast<std::size_t>(i)] * coords[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)];
      const Vec3 b = body_force_density(bl, rho[e], factor, xg);
      for (int i = 0; i < n; ++i)
        for (int d = 0; d < 3; ++d)
          f[static_cast<std::size_t>(nidx[static_cast<std::size_t>(i)]) * kDofsPerNode + static_cast<std::size_t>(d)] +=
              s.N[static_cast<std::size_t>(i)] * b[static_cast<std::size_t>(d)] * dv;
    }
  };

  if (elset == nullptr) {
    for (std::size_t e = 0; e < mesh.num_elements(); ++e) do_element(e);
  } else {
    for (const Index eid : *elset) {
      const Index ei = mesh.element_index(eid);
      if (ei >= 0) do_element(static_cast<std::size_t>(ei));
    }
  }
}

// Resolve a face's local node indices and kind across the element families.
// Tets: 4 triangular faces (T3/T6). Hexes: 6 quad faces (Q4/Q8). Wedges: F1/F2
// triangular end faces, F3/F4/F5 quad side faces. Midside nodes are included for
// the quadratic variants. (spec: heat-transfer-analysis / loads — surface loads.)
FaceKind face_kind_of(ElementType type, int face, std::vector<int>& fn) {
  const auto fill = [&](const int* row, int count) {
    fn.assign(static_cast<std::size_t>(count), 0);
    for (int i = 0; i < count; ++i) fn[static_cast<std::size_t>(i)] = row[i];
  };
  switch (type) {
    case ElementType::C3D4:
    case ElementType::C3D10: {
      if (face < 1 || face > 4) throw std::runtime_error("tet face must be 1..4");
      const int nf = (type == ElementType::C3D4) ? 3 : 6;
      fill(kFace[face - 1], nf);
      return {false, nf};
    }
    case ElementType::C3D8:
    case ElementType::C3D8R:
    case ElementType::C3D20:
    case ElementType::C3D20R: {
      if (face < 1 || face > 6) throw std::runtime_error("hex face must be 1..6");
      const bool lin = (type == ElementType::C3D8 || type == ElementType::C3D8R);
      const int nf = lin ? 4 : 8;
      fill(kFaceHex[face - 1], nf);
      return {true, nf};
    }
    case ElementType::C3D6:
    case ElementType::C3D15: {
      if (face < 1 || face > 5) throw std::runtime_error("wedge face must be 1..5");
      const bool tri = (face <= 2);
      const bool lin = (type == ElementType::C3D6);
      const int nf = tri ? (lin ? 3 : 6) : (lin ? 4 : 8);
      fill(kFaceWedge[face - 1], nf);
      return {!tri, nf};
    }
  }
  throw std::runtime_error("unsupported element type for surface load");
}

}  // namespace

std::vector<int> face_nodes(ElementType type, int face) {
  std::vector<int> fn;
  face_kind_of(type, face, fn);
  return fn;
}

FaceSurface face_surface_integrals(const Mesh& mesh, Index elem_index, int face) {
  const Element& elem = mesh.elements()[static_cast<std::size_t>(elem_index)];
  const FaceGeom fg = gather_face(mesh, elem, face);
  const int nf = fg.kind.nf;

  FaceSurface fs;
  fs.gnode = fg.gnode;
  fs.load.assign(static_cast<std::size_t>(nf), 0.0);
  fs.mass.assign(static_cast<std::size_t>(nf) * static_cast<std::size_t>(nf), 0.0);

  std::array<Real, 8> N{};
  Vec3 txi{}, teta{};
  for (const FaceGP& gp : face_gauss(fg.kind)) {
    face_point(fg, gp, N, txi, teta);
    const Vec3 nrm{txi[1] * teta[2] - txi[2] * teta[1],
                   txi[2] * teta[0] - txi[0] * teta[2],
                   txi[0] * teta[1] - txi[1] * teta[0]};
    const Real dA =
        std::sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]) * gp.w;
    for (int i = 0; i < nf; ++i) {
      const Real Ni = N[static_cast<std::size_t>(i)];
      fs.load[static_cast<std::size_t>(i)] += Ni * dA;
      for (int j = 0; j < nf; ++j)
        fs.mass[static_cast<std::size_t>(i) * static_cast<std::size_t>(nf) +
                static_cast<std::size_t>(j)] += Ni * N[static_cast<std::size_t>(j)] * dA;
    }
  }
  return fs;
}

// Emit one physical face point from a (xi,eta) sample with a natural-domain weight w.
FacePoint face_point_phys(const FaceGeom& fg, Real xi, Real eta, Real w) {
  std::array<Real, 8> N{};
  Vec3 txi{}, teta{};
  face_point(fg, FaceGP{xi, eta, 1.0}, N, txi, teta);
  const Vec3 nrm{txi[1] * teta[2] - txi[2] * teta[1],
                 txi[2] * teta[0] - txi[0] * teta[2],
                 txi[0] * teta[1] - txi[1] * teta[0]};
  const Real jac = std::sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]);
  Vec3 x{0, 0, 0};
  for (int i = 0; i < fg.kind.nf; ++i)
    for (int d = 0; d < 3; ++d)
      x[static_cast<std::size_t>(d)] += N[static_cast<std::size_t>(i)] *
                                        fg.x[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)];
  return FacePoint{x, normalize(nrm), jac * w};
}

std::vector<FacePoint> face_gauss_points(const Mesh& mesh, Index elem_index, int face,
                                         int sub) {
  const Element& elem = mesh.elements()[static_cast<std::size_t>(elem_index)];
  const FaceGeom fg = gather_face(mesh, elem, face);
  if (sub < 1) sub = 1;

  std::vector<FacePoint> pts;
  // Tile the face natural domain into sub x sub subcells and apply the base Gauss rule
  // in each; the weights scale by the subcell size so Σ w still equals the face area.
  // Quads live on [-1,1]^2 (subcell side 2/sub, weight factor 1/sub^2); triangles use
  // barycentric (xi,eta) on the unit triangle (uniform 1/sub^2 area split by an affine
  // map of the base rule into each of the sub^2 small triangles).
  const Real inv = 1.0 / static_cast<Real>(sub);
  if (fg.kind.quad) {
    for (int iu = 0; iu < sub; ++iu)
      for (int iv = 0; iv < sub; ++iv) {
        const Real u0 = -1.0 + 2.0 * iu * inv, v0 = -1.0 + 2.0 * iv * inv;
        for (const FaceGP& gp : face_gauss(fg.kind)) {
          const Real xi = u0 + (gp.xi + 1.0) * inv;   // map [-1,1] -> subcell
          const Real eta = v0 + (gp.eta + 1.0) * inv;
          pts.push_back(face_point_phys(fg, xi, eta, gp.w * inv * inv));
        }
      }
    return pts;
  }
  // Triangle: split the unit reference triangle into sub^2 small triangles (the usual
  // upright + inverted tiling) and map the base 3-point rule into each with weight/sub^2.
  for (int row = 0; row < sub; ++row)
    for (int col = 0; col < sub - row; ++col) {
      // Upright small triangle with corners (a0,b0),(a0+inv,b0),(a0,b0+inv).
      const Real a0 = col * inv, b0 = row * inv;
      for (const FaceGP& gp : face_gauss(fg.kind)) {
        const Real xi = a0 + gp.xi * inv, eta = b0 + gp.eta * inv;
        pts.push_back(face_point_phys(fg, xi, eta, gp.w * inv * inv));
      }
      if (col < sub - row - 1) {  // inverted small triangle fills the gap
        const Real ai = a0 + inv, bi = b0 + inv;
        for (const FaceGP& gp : face_gauss(fg.kind)) {
          const Real xi = ai - gp.xi * inv, eta = bi - gp.eta * inv;
          pts.push_back(face_point_phys(fg, xi, eta, gp.w * inv * inv));
        }
      }
    }
  return pts;
}

FaceFrame face_frame_at(const Mesh& mesh, Index elem_index, int face, Real xi,
                        Real eta) {
  const Element& elem = mesh.elements()[static_cast<std::size_t>(elem_index)];
  const FaceGeom fg = gather_face(mesh, elem, face);
  std::array<Real, 8> N{};
  Vec3 txi{}, teta{};
  face_point(fg, FaceGP{xi, eta, 1.0}, N, txi, teta);
  Vec3 x{0, 0, 0};
  for (int i = 0; i < fg.kind.nf; ++i)
    for (int d = 0; d < 3; ++d)
      x[static_cast<std::size_t>(d)] += N[static_cast<std::size_t>(i)] *
                                        fg.x[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)];
  return FaceFrame{x, txi, teta};
}

FaceEval face_eval_deformed(const Mesh& mesh, Index elem_index, int face, Real xi,
                            Real eta, const std::vector<Vec3>& u) {
  const Element& elem = mesh.elements()[static_cast<std::size_t>(elem_index)];
  FaceGeom fg = gather_face(mesh, elem, face);
  // Add the nodal displacement to each face node's reference coordinate.
  for (int i = 0; i < fg.kind.nf; ++i) {
    const std::size_t gi = static_cast<std::size_t>(fg.gnode[static_cast<std::size_t>(i)]);
    for (int d = 0; d < 3; ++d)
      fg.x[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)] +=
          u[gi][static_cast<std::size_t>(d)];
  }
  std::array<Real, 8> N{};
  Vec3 txi{}, teta{};
  face_point(fg, FaceGP{xi, eta, 1.0}, N, txi, teta);
  FaceEval ev;
  ev.gnode = fg.gnode;
  ev.N.assign(static_cast<std::size_t>(fg.kind.nf), 0.0);
  for (int i = 0; i < fg.kind.nf; ++i) {
    ev.N[static_cast<std::size_t>(i)] = N[static_cast<std::size_t>(i)];
    for (int d = 0; d < 3; ++d)
      ev.x[static_cast<std::size_t>(d)] += N[static_cast<std::size_t>(i)] *
                                           fg.x[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)];
  }
  ev.t_xi = txi;
  ev.t_eta = teta;
  return ev;
}

std::vector<Real> thermal_load_vector(const Model& model, Real lambda) {
  const Mesh& mesh = model.mesh;
  std::vector<Real> q(mesh.num_nodes(), 0.0);

  for (const Cflux& cf : model.cfluxes) {
    const Index ni = mesh.node_index(cf.node_id);
    if (ni >= 0) {
      const Real s = model.amplitude_factor(cf.amplitude, lambda);
      q[static_cast<std::size_t>(ni)] += s * cf.value;
    }
  }

  for (const Dflux& df : model.dfluxes) {
    const Index ei = mesh.element_index(df.elem_id);
    if (ei < 0) throw std::runtime_error("*DFLUX references unknown element");
    const Real s = model.amplitude_factor(df.amplitude, lambda);
    accumulate_flux_face(mesh, ei, df.face, s * df.flux, q);
  }

  return q;
}

std::vector<Real> thermal_load_vector(const Model& model) {
  return thermal_load_vector(model, 1.0);
}

std::vector<Real> thermal_strain_load_vector(const Model& model, Real lambda) {
  const Mesh& mesh = model.mesh;
  std::vector<Real> f(mesh.num_nodes() * kDofsPerNode, 0.0);
  if (model.applied_temperature.empty()) return f;

  const std::vector<ElasticIso> elastic = model.element_elastic();
  const std::vector<std::optional<Expansion>> expansion = model.element_expansion();
  const std::vector<bool> active = model.element_active_mask();
  std::vector<Vec3> coords;
  std::vector<Real> te;  // per-node temperature change (T - Tref), scaled by lambda

  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    if (!active[e]) continue;  // *MODEL CHANGE, REMOVE: no thermal-strain load
    const std::optional<Expansion>& exp = expansion[e];
    if (!exp || exp->empty()) continue;
    const Element& elem = mesh.elements()[e];
    const int n = nodes_per_element(elem.type);
    coords.resize(static_cast<std::size_t>(n));
    te.resize(static_cast<std::size_t>(n));
    bool any_temp = false;
    Real tsum = 0.0;  // sum of nodal applied temperatures (for alpha(T) evaluation)
    std::vector<Index> nidx(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
      const Index ni = mesh.node_index(elem.nodes[static_cast<std::size_t>(i)]);
      if (ni < 0) throw std::runtime_error("element references unknown node");
      nidx[static_cast<std::size_t>(i)] = ni;
      coords[static_cast<std::size_t>(i)] = mesh.nodes()[static_cast<std::size_t>(ni)].x;
      const auto it = model.applied_temperature.find(elem.nodes[static_cast<std::size_t>(i)]);
      const Real T = (it == model.applied_temperature.end()) ? exp->t_ref : it->second;
      te[static_cast<std::size_t>(i)] = lambda * (T - exp->t_ref);
      tsum += T;
      if (T != exp->t_ref) any_temp = true;
    }
    if (!any_temp) continue;
    // Temperature-dependent alpha(T) is evaluated at the element mean temperature
    // (constant tables return their single value, so a constant deck is unchanged).
    const Real alpha = exp->alpha.at(tsum / static_cast<Real>(n));
    const D6 D = elastic_iso_D(elastic[e]);
    const std::vector<Real> fe =
        element_thermal_load(elem.type, coords, D, alpha, te);
    for (int a = 0; a < n * kDofsPerNode; ++a)
      f[static_cast<std::size_t>(nidx[static_cast<std::size_t>(a / kDofsPerNode)]) *
            kDofsPerNode +
        static_cast<std::size_t>(a % kDofsPerNode)] += fe[static_cast<std::size_t>(a)];
  }
  return f;
}

std::vector<Real> thermal_strain_load_vector(const Model& model) {
  return thermal_strain_load_vector(model, 1.0);
}

std::vector<Real> external_load_vector(const Model& model, Real lambda) {
  const Mesh& mesh = model.mesh;
  std::vector<Real> f(mesh.num_nodes() * kDofsPerNode, 0.0);

  for (const Cload& cl : model.cloads) {
    const Index ni = mesh.node_index(cl.node_id);
    if (ni >= 0 && cl.comp >= 1 && cl.comp <= kDofsPerNode) {
      const Real s = model.amplitude_factor(cl.amplitude, lambda);
      f[static_cast<std::size_t>(ni) * kDofsPerNode + static_cast<std::size_t>(cl.comp - 1)] += s * cl.value;
    }
  }

  const std::vector<bool> active = model.element_active_mask();
  for (const Dload& dl : model.dloads) {
    const Index ei = mesh.element_index(dl.elem_id);
    if (ei < 0) throw std::runtime_error("*DLOAD references unknown element");
    // A pressure on a *MODEL CHANGE, REMOVE element carries no load (the face is gone).
    if (!active[static_cast<std::size_t>(ei)]) continue;
    const Real s = model.amplitude_factor(dl.amplitude, lambda);
    accumulate_pressure_face(mesh, ei, dl.face, s * dl.pressure, f);
  }

  for (const BodyLoad& bl : model.body_loads) {
    const Real s = model.amplitude_factor(bl.amplitude, lambda);
    accumulate_body_load(model, bl, s, f);
  }

  return f;
}

std::vector<Real> external_load_vector(const Model& model) {
  // Full magnitude: step fraction lambda = 1. Unamplified loads get factor 1 (the
  // linear-path behavior, independent of total time); amplitudes are sampled at
  // their end-of-step time.
  return external_load_vector(model, 1.0);
}

}  // namespace cxpp::fem
