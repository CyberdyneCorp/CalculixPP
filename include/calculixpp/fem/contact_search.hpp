#pragma once
#include <array>
#include <cstddef>
#include <utility>
#include <vector>

#include "calculixpp/core/mesh.hpp"
#include "calculixpp/core/types.hpp"

// Geometric contact-search engine (spec: contact-search — the NEW capability). This
// header is the geometric CORE only: spatial proximity search over master faces
// (task 1.1), node-to-surface closest-point projection with signed gap and a local
// normal/tangent frame (task 1.2), and a per-contact state struct that later carries
// pressure/area for mechanical AND thermal contact reuse (task 1.6). No system
// assembly here — the operator generation (task 1.5) that feeds the Newton tangent is
// a later workstream; this piece is unit-tested against analytic projections and does
// not touch any solve path, so a deck with no contact is byte-for-byte unchanged.
//
// Math reference (reimplemented, not copied): CalculiX gencontelem_n2f.f (the n2f
// projection) and the near2d spatial bucketing — algorithm reference only.
namespace cxpp::fem {

// A master face the search projects onto: an element face identified by the element's
// dense mesh index and its face number (1-based, the same numbering the *DFLUX/*FILM
// face loads use). Resolved from a *SURFACE (TYPE=ELEMENT) by expand_master_faces.
struct MasterFace {
  Index elem_index{-1};  // dense index into mesh.elements()
  int face{0};           // 1-based face number (tet 1..4, hex 1..6, wedge 1..5)
};

// Resolve the element faces of a named surface into MasterFaces, using the element
// INDICES (not ids) so downstream geometry lookups are O(1). A TYPE=ELEMENT surface
// contributes its (element id, face) pairs directly; a TYPE=NODE surface has no faces
// and yields an empty list (node surfaces are slave-side only in node-to-surface
// contact). Faces whose element id is not in the mesh are skipped. (spec: contact —
// *CONTACT PAIR consumes named surfaces; the search engine sees resolved faces.)
std::vector<MasterFace> expand_master_faces(const Mesh& mesh, const Surface& surface);

// Axis-aligned bounding box of one master face, INFLATED by `pad` on every side (the
// search/adjust distance). Used both to build the spatial grid and to prune the
// candidate set. `pad >= 0`.
struct Aabb {
  Vec3 lo{0, 0, 0};
  Vec3 hi{0, 0, 0};
  bool contains(const Vec3& p) const;
  bool overlaps(const Aabb& o) const;
};
Aabb face_aabb(const Mesh& mesh, const MasterFace& mf, Real pad = 0.0);

// Uniform-grid spatial index over the master faces, so a slave node finds its nearby
// candidate faces WITHOUT the O(n_slave * n_master) all-pairs test (spec: contact-search
// — large surfaces avoid all-pairs cost). Each face's padded AABB is rasterized into the
// cells it overlaps; a query gathers the faces in the cell(s) covering the query point
// plus the search radius. The cell size auto-sizes to the mean padded-face extent so the
// average occupancy stays O(1); a degenerate (empty) index answers every query with an
// empty candidate list.
class MasterGrid {
 public:
  // Build the grid over `faces` with a uniform search/adjust distance `search`. Faces
  // are stored by their index into `faces`; the caller keeps `faces` alive and maps a
  // returned index back to its MasterFace.
  MasterGrid(const Mesh& mesh, const std::vector<MasterFace>& faces, Real search);

  Real search_distance() const { return search_; }
  std::size_t num_faces() const { return num_faces_; }

  // Indices (into the `faces` vector passed to the constructor) of every master face
  // whose padded AABB could be within the search distance of `p`. A superset of the
  // truly-close faces (broad phase); the caller narrows it by projecting. Deterministic
  // order (ascending face index), de-duplicated.
  std::vector<int> candidates(const Vec3& p) const;

 private:
  const Mesh* mesh_;
  const std::vector<MasterFace>* faces_;
  Real search_{0.0};
  Real cell_{1.0};
  Vec3 origin_{0, 0, 0};
  std::array<int, 3> dims_{1, 1, 1};
  std::vector<std::vector<int>> cells_;  // cell -> face indices
  std::size_t num_faces_{0};

  int cell_of(const Vec3& p, std::array<int, 3>& ijk) const;
  int flat(int i, int j, int k) const;
};

// The result of projecting a slave node onto a master face (spec: contact-search —
// node-to-surface projection AND the shared state for thermal contact). It carries the
// closest point, the LOCAL face coordinates of that point, the signed normal gap, and
// the local normal/tangent frame; `pressure` and `area` are left at their defaults here
// (the geometric core does not evaluate the contact law — task 1.4 does) so the SAME
// struct is the interface geometry that mechanical contact and thermal gap conductance/
// heat generation both consume (task 1.6). `valid` is false when no candidate face
// admits a usable projection.
struct ContactState {
  bool valid{false};
  int master_face{-1};   // index into the faces vector; -1 when invalid
  Real xi{0.0}, eta{0.0}; // local (parametric) coords of the projection on the face
  Vec3 point{0, 0, 0};   // projected point on the master face (physical)
  Vec3 normal{0, 0, 0};  // unit outward master-face normal at the projection
  Vec3 tangent1{0, 0, 0};// unit surface tangent (t_xi normalized)
  Vec3 tangent2{0, 0, 0};// unit surface tangent completing the right-handed frame
  // Signed normal gap g = (x_slave - x_proj) . normal. g < 0 is penetration (the slave
  // node is on the material side of the master face); g > 0 is a positive clearance.
  Real gap{0.0};
  // Filled by the later contact-law / thermal workstreams; defaults keep the geometric
  // core law-agnostic. (task 1.6 — exposed for thermal-contact reuse.)
  Real pressure{0.0};
  Real area{0.0};
};

// Project one point onto ONE master face by the closest-point (Newton) iteration on the
// face natural coordinates, clamping to the face domain, and return the state (point,
// local coords, signed gap, normal/tangent frame). Always returns a state with
// valid = true (a face is always a well-defined surface); the gap sign tells penetration
// from clearance. (spec: contact-search — slave node projected to master face.)
ContactState project_onto_face(const Mesh& mesh, const MasterFace& mf,
                               const Vec3& slave);

// Project a slave point onto its NEAREST candidate master face found through the grid,
// within the grid's search distance. Runs the broad-phase grid query, projects onto each
// candidate, and keeps the projection with the smallest absolute gap (closest surface).
// Returns an invalid state (valid = false) when no candidate lies within the search
// distance. This is the per-slave-node entry point the contact operator (task 1.5) will
// call each increment. (spec: contact-search — spatial proximity + projection.)
ContactState project_nearest(const Mesh& mesh, const std::vector<MasterFace>& faces,
                             const MasterGrid& grid, const Vec3& slave);

}  // namespace cxpp::fem
