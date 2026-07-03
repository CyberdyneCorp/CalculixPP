#include "calculixpp/fem/contact_search.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "calculixpp/core/element.hpp"
#include "calculixpp/fem/loads.hpp"

namespace cxpp::fem {
namespace {

Vec3 sub(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Real dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
Real norm(const Vec3& a) { return std::sqrt(dot(a, a)); }
Vec3 normalized(const Vec3& a) {
  const Real n = norm(a);
  if (n <= 0.0) return {0, 0, 0};
  return {a[0] / n, a[1] / n, a[2] / n};
}

// True when the face is a triangular face (its parametric domain is the unit reference
// triangle l1 = 1-xi-eta >= 0); false for quad faces on [-1,1]^2. Tets: all faces
// triangular. Wedges: faces 1,2 triangular, 3,4,5 quad. Hexes: all quad.
bool is_triangle_face(ElementType type, int face) {
  switch (type) {
    case ElementType::C3D4:
    case ElementType::C3D10:
      return true;
    case ElementType::C3D6:
    case ElementType::C3D15:
      return face <= 2;
    default:
      return false;
  }
}

// Clamp parametric coords to the face domain: [-1,1]^2 for a quad, the unit triangle
// (xi>=0, eta>=0, xi+eta<=1) for a triangle. Keeps a projection that walks off the face
// pinned to the nearest boundary point in parameter space (a simple, robust clamp; the
// closest-point-on-a-general-quad edge case is a later refinement).
void clamp_to_face(bool tri, Real& xi, Real& eta) {
  if (!tri) {
    xi = std::clamp(xi, -1.0, 1.0);
    eta = std::clamp(eta, -1.0, 1.0);
    return;
  }
  xi = std::clamp(xi, 0.0, 1.0);
  eta = std::clamp(eta, 0.0, 1.0);
  const Real s = xi + eta;
  if (s > 1.0) {  // project onto the hypotenuse xi+eta=1
    const Real d = 0.5 * (s - 1.0);
    xi -= d;
    eta -= d;
    xi = std::clamp(xi, 0.0, 1.0);
    eta = std::clamp(eta, 0.0, 1.0);
  }
}

// A starting parametric point at the face centroid.
void face_center(bool tri, Real& xi, Real& eta) {
  if (tri) { xi = 1.0 / 3.0; eta = 1.0 / 3.0; }
  else { xi = 0.0; eta = 0.0; }
}

// Centroid of an element's nodes. The raw face tangent cross product t_xi x t_eta is
// NOT reliably outward across the element families (its sign follows the stored face
// node ordering, which the pressure kernel consumes with its own sign). For contact the
// master normal must point OUT of the master solid (toward the slave), so we orient it
// away from the element centroid — a topology-robust choice, the same reasoning the
// cavity-radiation patch uses.
Vec3 element_centroid(const Mesh& mesh, const Element& elem) {
  Vec3 c{0, 0, 0};
  for (const Index nid : elem.nodes) {
    const Index ni = mesh.node_index(nid);
    const Vec3& x = mesh.nodes()[static_cast<std::size_t>(ni)].x;
    for (int d = 0; d < 3; ++d) c[static_cast<std::size_t>(d)] += x[static_cast<std::size_t>(d)];
  }
  const Real inv = 1.0 / static_cast<Real>(elem.nodes.size());
  return {c[0] * inv, c[1] * inv, c[2] * inv};
}

}  // namespace

bool Aabb::contains(const Vec3& p) const {
  return p[0] >= lo[0] && p[0] <= hi[0] && p[1] >= lo[1] && p[1] <= hi[1] &&
         p[2] >= lo[2] && p[2] <= hi[2];
}

bool Aabb::overlaps(const Aabb& o) const {
  return lo[0] <= o.hi[0] && hi[0] >= o.lo[0] && lo[1] <= o.hi[1] &&
         hi[1] >= o.lo[1] && lo[2] <= o.hi[2] && hi[2] >= o.lo[2];
}

std::vector<MasterFace> expand_master_faces(const Mesh& mesh, const Surface& surface) {
  std::vector<MasterFace> out;
  if (surface.type != Surface::Type::Element) return out;
  out.reserve(surface.faces.size());
  for (const auto& [elem_id, face] : surface.faces) {
    const Index ei = mesh.element_index(elem_id);
    if (ei >= 0) out.push_back(MasterFace{ei, face});
  }
  return out;
}

Aabb face_aabb(const Mesh& mesh, const MasterFace& mf, Real pad) {
  const Element& elem = mesh.elements()[static_cast<std::size_t>(mf.elem_index)];
  const std::vector<int> fn = face_nodes(elem.type, mf.face);
  Aabb box;
  box.lo = {std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max(),
            std::numeric_limits<Real>::max()};
  box.hi = {std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::lowest(),
            std::numeric_limits<Real>::lowest()};
  for (const int ln : fn) {
    const Index ni = mesh.node_index(elem.nodes[static_cast<std::size_t>(ln)]);
    const Vec3& x = mesh.nodes()[static_cast<std::size_t>(ni)].x;
    for (int d = 0; d < 3; ++d) {
      box.lo[static_cast<std::size_t>(d)] =
          std::min(box.lo[static_cast<std::size_t>(d)], x[static_cast<std::size_t>(d)]);
      box.hi[static_cast<std::size_t>(d)] =
          std::max(box.hi[static_cast<std::size_t>(d)], x[static_cast<std::size_t>(d)]);
    }
  }
  for (int d = 0; d < 3; ++d) {
    box.lo[static_cast<std::size_t>(d)] -= pad;
    box.hi[static_cast<std::size_t>(d)] += pad;
  }
  return box;
}

MasterGrid::MasterGrid(const Mesh& mesh, const std::vector<MasterFace>& faces,
                       Real search)
    : mesh_(&mesh), faces_(&faces), search_(std::max(search, 0.0)),
      num_faces_(faces.size()) {
  if (faces.empty()) return;

  // Overall padded extent and the mean face size -> a cell size that keeps ~O(1) faces
  // per cell. The cell must be at least the search distance so a single ring of cells
  // around the query point already covers the search radius.
  Aabb bounds = face_aabb(mesh, faces[0], search_);
  Real mean_extent = 0.0;
  for (const MasterFace& mf : faces) {
    const Aabb b = face_aabb(mesh, mf, 0.0);
    bounds.lo[0] = std::min(bounds.lo[0], b.lo[0] - search_);
    bounds.lo[1] = std::min(bounds.lo[1], b.lo[1] - search_);
    bounds.lo[2] = std::min(bounds.lo[2], b.lo[2] - search_);
    bounds.hi[0] = std::max(bounds.hi[0], b.hi[0] + search_);
    bounds.hi[1] = std::max(bounds.hi[1], b.hi[1] + search_);
    bounds.hi[2] = std::max(bounds.hi[2], b.hi[2] + search_);
    mean_extent += (b.hi[0] - b.lo[0]) + (b.hi[1] - b.lo[1]) + (b.hi[2] - b.lo[2]);
  }
  mean_extent /= static_cast<Real>(3 * faces.size());
  cell_ = std::max({mean_extent, search_, 1e-12});
  origin_ = bounds.lo;
  for (int d = 0; d < 3; ++d) {
    const Real span = bounds.hi[static_cast<std::size_t>(d)] -
                      bounds.lo[static_cast<std::size_t>(d)];
    dims_[static_cast<std::size_t>(d)] =
        std::max(1, static_cast<int>(std::floor(span / cell_)) + 1);
  }
  cells_.assign(static_cast<std::size_t>(dims_[0]) * dims_[1] * dims_[2], {});

  // Rasterize each padded face AABB into the cells it overlaps.
  for (int f = 0; f < static_cast<int>(faces.size()); ++f) {
    const Aabb b = face_aabb(mesh, faces[static_cast<std::size_t>(f)], search_);
    std::array<int, 3> lo{}, hi{};
    for (int d = 0; d < 3; ++d) {
      lo[static_cast<std::size_t>(d)] = std::clamp(
          static_cast<int>(std::floor((b.lo[static_cast<std::size_t>(d)] -
                                       origin_[static_cast<std::size_t>(d)]) / cell_)),
          0, dims_[static_cast<std::size_t>(d)] - 1);
      hi[static_cast<std::size_t>(d)] = std::clamp(
          static_cast<int>(std::floor((b.hi[static_cast<std::size_t>(d)] -
                                       origin_[static_cast<std::size_t>(d)]) / cell_)),
          0, dims_[static_cast<std::size_t>(d)] - 1);
    }
    for (int i = lo[0]; i <= hi[0]; ++i)
      for (int j = lo[1]; j <= hi[1]; ++j)
        for (int k = lo[2]; k <= hi[2]; ++k)
          cells_[static_cast<std::size_t>(flat(i, j, k))].push_back(f);
  }
}

int MasterGrid::flat(int i, int j, int k) const {
  return (i * dims_[1] + j) * dims_[2] + k;
}

int MasterGrid::cell_of(const Vec3& p, std::array<int, 3>& ijk) const {
  for (int d = 0; d < 3; ++d)
    ijk[static_cast<std::size_t>(d)] = static_cast<int>(std::floor(
        (p[static_cast<std::size_t>(d)] - origin_[static_cast<std::size_t>(d)]) / cell_));
  return 0;
}

std::vector<int> MasterGrid::candidates(const Vec3& p) const {
  std::vector<int> out;
  if (cells_.empty()) return out;
  std::array<int, 3> c{};
  cell_of(p, c);
  // The query point plus its search radius spans at most one cell in each direction
  // (cell_ >= search_), so a 3x3x3 neighborhood of cells covers every face whose padded
  // AABB could reach the point.
  for (int di = -1; di <= 1; ++di)
    for (int dj = -1; dj <= 1; ++dj)
      for (int dk = -1; dk <= 1; ++dk) {
        const int i = c[0] + di, j = c[1] + dj, k = c[2] + dk;
        if (i < 0 || i >= dims_[0] || j < 0 || j >= dims_[1] || k < 0 || k >= dims_[2])
          continue;
        for (const int f : cells_[static_cast<std::size_t>(flat(i, j, k))])
          out.push_back(f);
      }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

ContactState project_onto_face(const Mesh& mesh, const MasterFace& mf,
                               const Vec3& slave) {
  const Element& elem = mesh.elements()[static_cast<std::size_t>(mf.elem_index)];
  const bool tri = is_triangle_face(elem.type, mf.face);

  Real xi = 0.0, eta = 0.0;
  face_center(tri, xi, eta);

  // Gauss-Newton closest-point iteration: solve r_a = (slave - x(xi,eta)) . t_a = 0 for
  // a = xi,eta, dropping the surface-curvature (second-derivative) term of the Jacobian
  // so J_ab = -t_a . t_b (the standard robust contact projection). Clamp to the face
  // domain after each step so the point stays on the face.
  FaceFrame fr = face_frame_at(mesh, mf.elem_index, mf.face, xi, eta);
  for (int it = 0; it < 20; ++it) {
    fr = face_frame_at(mesh, mf.elem_index, mf.face, xi, eta);
    const Vec3 d = sub(slave, fr.x);
    const Real r0 = dot(d, fr.t_xi), r1 = dot(d, fr.t_eta);
    const Real a = dot(fr.t_xi, fr.t_xi), b = dot(fr.t_xi, fr.t_eta),
               cc = dot(fr.t_eta, fr.t_eta);
    const Real det = a * cc - b * b;
    if (std::fabs(det) < 1e-30) break;  // degenerate tangent basis
    // Solve [a b; b cc] [dxi; deta] = [r0; r1].
    const Real dxi = (cc * r0 - b * r1) / det;
    const Real deta = (a * r1 - b * r0) / det;
    xi += dxi;
    eta += deta;
    clamp_to_face(tri, xi, eta);
    if (std::fabs(dxi) + std::fabs(deta) < 1e-13) break;
  }
  fr = face_frame_at(mesh, mf.elem_index, mf.face, xi, eta);

  ContactState st;
  st.valid = true;
  st.xi = xi;
  st.eta = eta;
  st.point = fr.x;
  Vec3 nrm = normalized(cross(fr.t_xi, fr.t_eta));
  // Orient outward: the master normal points from the element centroid toward the face.
  const Vec3 outward = sub(fr.x, element_centroid(mesh, elem));
  if (dot(nrm, outward) < 0.0) nrm = {-nrm[0], -nrm[1], -nrm[2]};
  st.normal = nrm;
  st.tangent1 = normalized(fr.t_xi);
  st.tangent2 = normalized(cross(st.normal, st.tangent1));
  st.gap = dot(sub(slave, fr.x), st.normal);
  return st;
}

ContactState project_nearest(const Mesh& mesh, const std::vector<MasterFace>& faces,
                             const MasterGrid& grid, const Vec3& slave) {
  ContactState best;  // valid=false by default
  Real best_dist = std::numeric_limits<Real>::max();
  for (const int f : grid.candidates(slave)) {
    ContactState st = project_onto_face(mesh, faces[static_cast<std::size_t>(f)], slave);
    const Real dist = norm(sub(slave, st.point));  // true 3D distance to the surface
    if (dist > grid.search_distance()) continue;   // outside the search/adjust distance
    if (dist < best_dist) {
      best_dist = dist;
      st.master_face = f;
      best = st;
    }
  }
  return best;
}

}  // namespace cxpp::fem
