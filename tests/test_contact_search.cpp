// Contact-search geometric core (spec: contact-search — tasks 1.1, 1.2, 1.6). These
// tests validate the geometry ANALYTICALLY (no solve path, no deck behavior): a slave
// node projected onto a KNOWN planar master face lands at the right point with the
// right signed gap; a node exactly on the face has g ~ 0; the spatial proximity search
// returns only nearby faces and scales sub-quadratically (it never visits every master
// face). The per-contact ContactState struct exposes gap/pressure/area/projection for
// the later thermal-contact reuse (task 1.6).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

#include "calculixpp/core/mesh.hpp"
#include "calculixpp/fem/contact_search.hpp"
#include "check.hpp"

using namespace cxpp;
using namespace cxpp::fem;

namespace {

Real dist(const Vec3& a, const Vec3& b) {
  const Real dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
  return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// A unit C3D8 hex on [0,1]^3. CalculiX C3D8 node order: bottom face (z=0) CCW 1-4,
// top face (z=1) CCW 5-8. Its face F2 (nodes 5,8,7,6) is the top z=1 plane with
// outward normal +z. Returns a mesh containing exactly this element (element id 1).
Mesh unit_hex() {
  Mesh mesh;
  mesh.add_node(1, {0, 0, 0});
  mesh.add_node(2, {1, 0, 0});
  mesh.add_node(3, {1, 1, 0});
  mesh.add_node(4, {0, 1, 0});
  mesh.add_node(5, {0, 0, 1});
  mesh.add_node(6, {1, 0, 1});
  mesh.add_node(7, {1, 1, 1});
  mesh.add_node(8, {0, 1, 1});
  mesh.add_element(1, ElementType::C3D8, {1, 2, 3, 4, 5, 6, 7, 8});
  return mesh;
}

// A grid of unit hexes tiled along +x: element e occupies x in [e-1, e]. Each element
// exposes its top face F2 (z=1, +z normal). Returns the mesh; `n` elements, ids 1..n.
Mesh hex_strip(int n) {
  Mesh mesh;
  for (int e = 0; e < n; ++e) {
    const Real x0 = static_cast<Real>(e);
    const Index b = static_cast<Index>(8 * e);
    mesh.add_node(b + 1, {x0, 0, 0});
    mesh.add_node(b + 2, {x0 + 1, 0, 0});
    mesh.add_node(b + 3, {x0 + 1, 1, 0});
    mesh.add_node(b + 4, {x0, 1, 0});
    mesh.add_node(b + 5, {x0, 0, 1});
    mesh.add_node(b + 6, {x0 + 1, 0, 1});
    mesh.add_node(b + 7, {x0 + 1, 1, 1});
    mesh.add_node(b + 8, {x0, 1, 1});
    mesh.add_element(static_cast<Index>(e + 1), ElementType::C3D8,
                     {b + 1, b + 2, b + 3, b + 4, b + 5, b + 6, b + 7, b + 8});
  }
  return mesh;
}

// (1.2) Project a node above the top face onto the z=1 plane: closest point is the
// (x,y) footprint at z=1, the outward normal is +z, and a node ABOVE the face has a
// positive gap (clearance) while a node BELOW (inside the block) has a negative gap
// (penetration).
void test_projection_planar_face() {
  const Mesh mesh = unit_hex();
  const MasterFace top{mesh.element_index(1), 2};  // z=1 face, outward +z

  const Vec3 above{0.3, 0.7, 1.25};  // 0.25 above the face
  ContactState s = project_onto_face(mesh, top, above);
  CX_CHECK(s.valid);
  CX_NEAR(s.point[0], 0.3, 1e-10);
  CX_NEAR(s.point[1], 0.7, 1e-10);
  CX_NEAR(s.point[2], 1.0, 1e-10);
  CX_NEAR(s.normal[0], 0.0, 1e-10);
  CX_NEAR(s.normal[1], 0.0, 1e-10);
  CX_NEAR(s.normal[2], 1.0, 1e-10);
  CX_NEAR(s.gap, 0.25, 1e-10);  // above -> positive clearance

  // Orthonormal right-handed frame.
  CX_NEAR(s.tangent1[0] * s.normal[0] + s.tangent1[1] * s.normal[1] +
              s.tangent1[2] * s.normal[2], 0.0, 1e-10);
  CX_NEAR(s.tangent2[0] * s.normal[0] + s.tangent2[1] * s.normal[1] +
              s.tangent2[2] * s.normal[2], 0.0, 1e-10);
  CX_NEAR(dist(s.tangent1, {0, 0, 0}), 1.0, 1e-10);
  CX_NEAR(dist(s.tangent2, {0, 0, 0}), 1.0, 1e-10);

  const Vec3 below{0.3, 0.7, 0.9};  // 0.1 into the block -> penetration
  ContactState p = project_onto_face(mesh, top, below);
  CX_NEAR(p.point[2], 1.0, 1e-10);
  CX_NEAR(p.gap, -0.1, 1e-10);  // below the face -> negative gap (penetration)
  CX_CHECK(p.gap < 0.0);
}

// (1.2) A node exactly on the face has g ~ 0 and projects to itself.
void test_node_on_face_zero_gap() {
  const Mesh mesh = unit_hex();
  const MasterFace top{mesh.element_index(1), 2};
  const Vec3 on{0.5, 0.5, 1.0};
  ContactState s = project_onto_face(mesh, top, on);
  CX_NEAR(s.gap, 0.0, 1e-10);
  CX_NEAR(dist(s.point, on), 0.0, 1e-10);
}

// (1.2) Projection clamps to the face: a slave node whose footprint is OUTSIDE the face
// footprint projects to the nearest point on the face boundary, not off it.
void test_projection_clamps_to_face() {
  const Mesh mesh = unit_hex();
  const MasterFace top{mesh.element_index(1), 2};
  const Vec3 outside{2.0, 0.5, 1.1};  // x=2 is beyond the [0,1] face
  ContactState s = project_onto_face(mesh, top, outside);
  CX_NEAR(s.point[0], 1.0, 1e-9);  // clamped to the x=1 edge of the face
  CX_NEAR(s.point[1], 0.5, 1e-9);
  CX_NEAR(s.point[2], 1.0, 1e-9);
}

// (1.1) Proximity search returns ONLY nearby faces. A strip of 40 hexes, each exposing
// its top face; a slave node near element 20 finds candidates only from the immediate
// neighborhood (a bounded constant), never all 40 faces, and the nearest projection is
// onto element 20's top face.
void test_proximity_returns_only_nearby() {
  const int n = 40;
  const Mesh mesh = hex_strip(n);
  std::vector<MasterFace> faces;
  for (int e = 1; e <= n; ++e) faces.push_back({mesh.element_index(e), 2});

  const Real search = 0.5;
  MasterGrid grid(mesh, faces, search);
  CX_CHECK(grid.num_faces() == static_cast<std::size_t>(n));

  // A node just above the center of element 20 (x in [19,20]).
  const Vec3 slave{19.5, 0.5, 1.1};
  std::vector<int> cand = grid.candidates(slave);
  CX_CHECK(!cand.empty());
  CX_CHECK(cand.size() < static_cast<std::size_t>(n));  // NOT all faces (no all-pairs)
  // Every candidate face's padded AABB must lie near the query in x (the strip axis):
  // its x-range overlaps [slave_x - cell, slave_x + cell]. Far faces are pruned.
  for (const int f : cand) {
    const Aabb b = face_aabb(mesh, faces[static_cast<std::size_t>(f)], search);
    CX_CHECK(b.hi[0] >= slave[0] - 3.0 && b.lo[0] <= slave[0] + 3.0);
  }
  // A far-away face (element 1, x in [0,1]) must NOT be a candidate.
  const int far = 0;  // faces[0] is element 1
  CX_CHECK(std::find(cand.begin(), cand.end(), far) == cand.end());

  // The nearest projection is onto element 20's top face at (19.5, 0.5, 1).
  ContactState s = project_nearest(mesh, faces, grid, slave);
  CX_CHECK(s.valid);
  CX_CHECK(s.master_face == 19);  // faces[19] == element 20
  CX_NEAR(s.point[0], 19.5, 1e-9);
  CX_NEAR(s.point[1], 0.5, 1e-9);
  CX_NEAR(s.point[2], 1.0, 1e-9);
  CX_NEAR(s.gap, 0.1, 1e-9);
}

// (1.1) Scaling: the total candidate work across every slave node stays LINEAR in the
// surface size, not quadratic. For a strip of n master faces we query one slave node
// over each face and assert the summed candidate count is O(n) (bounded by a small
// constant per query), proving the grid avoids the O(n^2) all-pairs comparison.
void test_proximity_scales() {
  for (const int n : {20, 40, 80}) {
    const Mesh mesh = hex_strip(n);
    std::vector<MasterFace> faces;
    for (int e = 1; e <= n; ++e) faces.push_back({mesh.element_index(e), 2});
    MasterGrid grid(mesh, faces, 0.5);

    std::size_t total = 0;
    for (int e = 0; e < n; ++e) {
      const Vec3 slave{static_cast<Real>(e) + 0.5, 0.5, 1.05};
      total += grid.candidates(slave).size();
    }
    // Each query touches a bounded neighborhood, so total <= C*n with a small C. A true
    // all-pairs scan would be n*n; assert we are far below that.
    CX_CHECK(total <= static_cast<std::size_t>(12 * n));
    CX_CHECK(total < static_cast<std::size_t>(n) * static_cast<std::size_t>(n));
  }
}

// (1.1) A query outside the whole surface (beyond the search distance) yields no
// candidates and an invalid nearest projection.
void test_proximity_out_of_range() {
  const Mesh mesh = hex_strip(5);
  std::vector<MasterFace> faces;
  for (int e = 1; e <= 5; ++e) faces.push_back({mesh.element_index(e), 2});
  MasterGrid grid(mesh, faces, 0.2);

  const Vec3 far{2.5, 0.5, 5.0};  // 4 units above every face, search only 0.2
  CX_CHECK(grid.candidates(far).empty());
  ContactState s = project_nearest(mesh, faces, grid, far);
  CX_CHECK(!s.valid);
}

// (1.6) The ContactState carries the projection geometry AND placeholder pressure/area
// so thermal contact can reuse the same struct. The geometric core leaves them zero.
void test_state_exposes_thermal_fields() {
  const Mesh mesh = unit_hex();
  const MasterFace top{mesh.element_index(1), 2};
  ContactState s = project_onto_face(mesh, top, {0.5, 0.5, 1.2});
  CX_NEAR(s.pressure, 0.0, 0.0);  // law-agnostic geometric core
  CX_NEAR(s.area, 0.0, 0.0);
  CX_CHECK(s.master_face == -1 || s.master_face >= 0);  // set by project_nearest
}

// expand_master_faces resolves a *SURFACE (element type) into MasterFaces by element
// index, and yields nothing for a node surface.
void test_expand_master_faces() {
  const Mesh mesh = hex_strip(3);
  Surface surf;
  surf.name = "MASTER";
  surf.type = Surface::Type::Element;
  surf.faces = {{1, 2}, {2, 2}, {3, 2}};
  std::vector<MasterFace> mf = expand_master_faces(mesh, surf);
  CX_CHECK(mf.size() == 3);
  CX_CHECK(mf[0].elem_index == mesh.element_index(1));
  CX_CHECK(mf[0].face == 2);

  Surface nsurf;
  nsurf.type = Surface::Type::Node;
  nsurf.nodes = {1, 2, 3};
  CX_CHECK(expand_master_faces(mesh, nsurf).empty());
}

}  // namespace

int main() {
  test_projection_planar_face();
  test_node_on_face_zero_gap();
  test_projection_clamps_to_face();
  test_proximity_returns_only_nearby();
  test_proximity_scales();
  test_proximity_out_of_range();
  test_state_exposes_thermal_fields();
  test_expand_master_faces();
  CX_MAIN_RETURN();
}
