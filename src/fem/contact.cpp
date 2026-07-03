#include "calculixpp/fem/contact.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>

#include "calculixpp/core/element.hpp"
#include "calculixpp/fem/constraint_transform.hpp"
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

void clamp_to_face(bool tri, Real& xi, Real& eta) {
  if (!tri) {
    xi = std::clamp(xi, -1.0, 1.0);
    eta = std::clamp(eta, -1.0, 1.0);
    return;
  }
  xi = std::clamp(xi, 0.0, 1.0);
  eta = std::clamp(eta, 0.0, 1.0);
  const Real s = xi + eta;
  if (s > 1.0) {
    const Real d = 0.5 * (s - 1.0);
    xi -= d;
    eta -= d;
    xi = std::clamp(xi, 0.0, 1.0);
    eta = std::clamp(eta, 0.0, 1.0);
  }
}

// Centroid of an element's nodes in the DEFORMED configuration (reference + u), used to
// orient the master normal outward (toward the slave) — the same topology-robust choice
// the geometric search uses on the reference mesh.
Vec3 deformed_centroid(const Mesh& mesh, const Element& elem,
                       const std::vector<Vec3>& u) {
  Vec3 c{0, 0, 0};
  for (const Index nid : elem.nodes) {
    const std::size_t ni = static_cast<std::size_t>(mesh.node_index(nid));
    const Vec3& x = mesh.nodes()[ni].x;
    for (int d = 0; d < 3; ++d)
      c[static_cast<std::size_t>(d)] += x[static_cast<std::size_t>(d)] + u[ni][static_cast<std::size_t>(d)];
  }
  const Real inv = 1.0 / static_cast<Real>(elem.nodes.size());
  return {c[0] * inv, c[1] * inv, c[2] * inv};
}

// The result of projecting a deformed slave point onto a deformed master face.
struct Projection {
  FaceEval ev;        // deformed face node ids + shape values + frame at the projection
  Vec3 normal{0, 0, 0};
  Real gap{0.0};      // (x_slave - x_proj) . normal; g < 0 = penetration
  Real dist{0.0};     // |x_slave - x_proj|, for nearest-face selection
};

// Closest-point projection of `slave` onto master face `mf` in the deformed config,
// mirroring project_onto_face but on (reference + u) geometry. Returns the deformed face
// evaluation at the footprint, the outward normal, and the signed gap.
Projection project_deformed(const Mesh& mesh, const MasterFace& mf, const Vec3& slave,
                            const std::vector<Vec3>& u) {
  const Element& elem = mesh.elements()[static_cast<std::size_t>(mf.elem_index)];
  const bool tri = is_triangle_face(elem.type, mf.face);
  Real xi = tri ? 1.0 / 3.0 : 0.0;
  Real eta = tri ? 1.0 / 3.0 : 0.0;

  FaceEval ev = face_eval_deformed(mesh, mf.elem_index, mf.face, xi, eta, u);
  for (int it = 0; it < 20; ++it) {
    ev = face_eval_deformed(mesh, mf.elem_index, mf.face, xi, eta, u);
    const Vec3 d = sub(slave, ev.x);
    const Real r0 = dot(d, ev.t_xi), r1 = dot(d, ev.t_eta);
    const Real a = dot(ev.t_xi, ev.t_xi), b = dot(ev.t_xi, ev.t_eta),
               cc = dot(ev.t_eta, ev.t_eta);
    const Real det = a * cc - b * b;
    if (std::fabs(det) < 1e-30) break;
    const Real dxi = (cc * r0 - b * r1) / det;
    const Real deta = (a * r1 - b * r0) / det;
    xi += dxi;
    eta += deta;
    clamp_to_face(tri, xi, eta);
    if (std::fabs(dxi) + std::fabs(deta) < 1e-13) break;
  }
  ev = face_eval_deformed(mesh, mf.elem_index, mf.face, xi, eta, u);

  Projection pr;
  pr.ev = ev;
  Vec3 nrm = normalized(cross(ev.t_xi, ev.t_eta));
  const Vec3 outward = sub(ev.x, deformed_centroid(mesh, elem, u));
  if (dot(nrm, outward) < 0.0) nrm = {-nrm[0], -nrm[1], -nrm[2]};
  pr.normal = nrm;
  pr.gap = dot(sub(slave, ev.x), nrm);
  pr.dist = norm(sub(slave, ev.x));
  return pr;
}

// Representative element size (max edge-ish extent) of the master element, for the auto
// hard-penalty estimate kappa ~ c * E * h.
Real element_size(const Mesh& mesh, const Element& elem) {
  Vec3 lo{std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max(),
          std::numeric_limits<Real>::max()};
  Vec3 hi{std::numeric_limits<Real>::lowest(), std::numeric_limits<Real>::lowest(),
          std::numeric_limits<Real>::lowest()};
  for (const Index nid : elem.nodes) {
    const Vec3& x = mesh.nodes()[static_cast<std::size_t>(mesh.node_index(nid))].x;
    for (int d = 0; d < 3; ++d) {
      lo[static_cast<std::size_t>(d)] = std::min(lo[static_cast<std::size_t>(d)], x[static_cast<std::size_t>(d)]);
      hi[static_cast<std::size_t>(d)] = std::max(hi[static_cast<std::size_t>(d)], x[static_cast<std::size_t>(d)]);
    }
  }
  return std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2], 1e-12});
}

// Resolve the slave nodes of a *CONTACT PAIR's slave surface: a TYPE=NODE surface gives
// its node list; a TYPE=ELEMENT surface contributes the (de-duplicated) nodes of its
// faces. Returned as dense mesh node indices.
std::vector<Index> resolve_slave_nodes(const Mesh& mesh, const Surface& surf) {
  std::vector<Index> out;
  if (surf.type == Surface::Type::Node) {
    for (const Index nid : surf.nodes) {
      const Index ni = mesh.node_index(nid);
      if (ni >= 0) out.push_back(ni);
    }
  } else {
    for (const auto& [elem_id, face] : surf.faces) {
      const Index ei = mesh.element_index(elem_id);
      if (ei < 0) continue;
      const Element& elem = mesh.elements()[static_cast<std::size_t>(ei)];
      for (const int ln : face_nodes(elem.type, face)) {
        const Index ni = mesh.node_index(elem.nodes[static_cast<std::size_t>(ln)]);
        if (ni >= 0) out.push_back(ni);
      }
    }
  }
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
  return out;
}

// Effective penalty stiffness (force per unit penetration) for a resolved pair. Linear
// uses the user slope; hard with an explicit k uses it; hard with k == 0 estimates
// kappa ~ 50 * E * h from the first master element (E from its material), a large
// spring giving near-zero overclosure. Exponential carries its own parameters and this
// returns p0/c0 as a fallback slope near contact. (spec: contact — surface behavior.)
Real resolve_penalty(const Model& model, const SurfaceBehavior& beh,
                     const std::vector<MasterFace>& faces) {
  if (beh.law == SurfaceBehavior::Law::Linear || beh.law == SurfaceBehavior::Law::Hard) {
    if (beh.k > 0.0) return beh.k;
  }
  if (beh.law == SurfaceBehavior::Law::Exponential && beh.c0 > 0.0)
    return beh.p0 / beh.c0;  // slope of the softened law near g = 0
  // Auto hard estimate from the first master element's stiffness and size.
  Real e_mod = 210000.0, h = 1.0;
  if (!faces.empty()) {
    const Element& elem =
        model.mesh.elements()[static_cast<std::size_t>(faces.front().elem_index)];
    h = element_size(model.mesh, elem);
    const std::vector<ElasticIso> ei = model.element_elastic();
    const std::size_t idx = static_cast<std::size_t>(faces.front().elem_index);
    if (idx < ei.size() && ei[idx].E > 0.0) e_mod = ei[idx].E;
  }
  return 50.0 * e_mod * h;
}

}  // namespace

Real contact_force_weight(const SurfaceBehavior& beh, Real kappa, Real g) {
  if (g >= 0.0) return 0.0;  // clearance -> no contact
  const Real pen = -g;       // penetration depth (> 0)
  switch (beh.law) {
    case SurfaceBehavior::Law::Hard:
    case SurfaceBehavior::Law::Linear:
      return kappa * pen;  // linear spring: F = kappa * penetration
    case SurfaceBehavior::Law::Exponential:
      // Softened exponential overclosure: force grows with penetration but saturates the
      // slope. F = kappa * c0 * (exp(pen/c0) - 1) reduces to kappa*pen for small pen and
      // stiffens deeper, a bounded-slope-at-zero variant of the CalculiX exponential law.
      if (beh.c0 > 0.0) return kappa * beh.c0 * (std::exp(pen / beh.c0) - 1.0);
      return kappa * pen;
  }
  return 0.0;
}

Vec3 friction_return_map(const Vec3& t_trial, Real mu, Real p_normal, bool& slipping) {
  slipping = false;
  const Real limit = mu * p_normal;  // Coulomb cap on tangential traction magnitude
  const Real mag = std::sqrt(dot(t_trial, t_trial));
  if (limit <= 0.0) {  // no normal pressure or mu == 0 -> no tangential traction can be held
    slipping = mag > 0.0;
    return {0, 0, 0};
  }
  if (mag <= limit) return t_trial;  // inside the cone -> stick, traction unchanged
  slipping = true;                   // on the cone -> slip, scale down to the limit
  const Real s = limit / mag;
  return {t_trial[0] * s, t_trial[1] * s, t_trial[2] * s};
}

std::vector<ResolvedContactPair> build_contact_pairs(const Model& model) {
  std::vector<ResolvedContactPair> out;
  for (const ContactPair& pair : model.contact_pairs) {
    if (pair.formulation != ContactPair::Formulation::NodeToSurface)
      throw std::runtime_error(
          "surface-to-surface (mortar) contact is not implemented yet; use "
          "TYPE=NODE TO SURFACE (tasks.md 1.3/2.1 mortar deferred)");
    const Surface* slave = model.mesh.surface(pair.slave_surface);
    const Surface* master = model.mesh.surface(pair.master_surface);
    if (slave == nullptr)
      throw std::runtime_error("*CONTACT PAIR references unknown slave surface '" +
                               pair.slave_surface + "'");
    if (master == nullptr)
      throw std::runtime_error("*CONTACT PAIR references unknown master surface '" +
                               pair.master_surface + "'");
    std::vector<MasterFace> faces = expand_master_faces(model.mesh, *master);
    if (faces.empty())
      throw std::runtime_error("*CONTACT PAIR master surface '" + pair.master_surface +
                               "' has no element faces (master must be TYPE=ELEMENT)");
    const SurfaceBehavior beh = model.contact_behavior(pair);
    const Friction fric = model.contact_friction(pair);
    // Default search distance: a few master-element sizes so a slave a bit off the face
    // still finds it (small-sliding). Explicit ADJUST/search overrides.
    Real search = pair.search;
    if (search <= 0.0)
      search = 2.0 * element_size(
                         model.mesh,
                         model.mesh.elements()[static_cast<std::size_t>(faces.front().elem_index)]);
    ResolvedContactPair rp{resolve_slave_nodes(model.mesh, *slave), std::move(faces),
                           MasterGrid(model.mesh, {}, search), beh, fric};
    // Rebuild the grid over the resolved faces (moved above) and finalize the penalty.
    rp.grid = MasterGrid(model.mesh, rp.master_faces, search);
    rp.penalty = resolve_penalty(model, beh, rp.master_faces);
    rp.clearance = pair.has_clearance ? pair.clearance : 0.0;
    // Tangential stick stiffness: the *FRICTION stick (lambda) value when given. Otherwise
    // an auto default that keeps the stick spring stiff enough to hold in the stick regime
    // yet moderate enough that the stick<->slip transition converges (a very large stick
    // stiffness relative to the material makes the transition increment chatter). We tie
    // the auto value to the master element's shear-scale stiffness (~ E*h capped by the
    // normal penalty), the same order as the elastic path the friction couples to.
    rp.stick = fric.stick;
    if (rp.stick <= 0.0) {
      Real e_mod = 210000.0, h = 1.0;
      const Element& melem =
          model.mesh.elements()[static_cast<std::size_t>(rp.master_faces.front().elem_index)];
      h = element_size(model.mesh, melem);
      const std::vector<ElasticIso> ei = model.element_elastic();
      const std::size_t idx = static_cast<std::size_t>(rp.master_faces.front().elem_index);
      if (idx < ei.size() && ei[idx].E > 0.0) e_mod = ei[idx].E;
      rp.stick = std::min(rp.penalty, e_mod * h * 0.1);  // ~ shear-scale, transition-stable
    }
    // Near-contact stiffness band: a small fraction of the master element size, so a
    // slave at first touch (gap ~ 0) already carries the penalty stiffness and the
    // interface has no rigid-body mode. Small enough not to glue truly-separated bodies.
    rp.adjust = 1e-3 * element_size(
                           model.mesh,
                           model.mesh.elements()[static_cast<std::size_t>(
                               rp.master_faces.front().elem_index)]);
    out.push_back(std::move(rp));
  }
  return out;
}

namespace {

// Scatter a penalty block into the reduced COO system through the constraint transform.
// `gdofs` are the participating global DOFs (slave x,y,z then each master node x,y,z);
// `K` is the dense m x m block (row-major). Mirrors fem::scatter_block but writes
// straight into sys.rows/cols/vals (the tangent is already flattened by the time the
// driver adds contact) and moves prescribed-DOF columns to sys.rhs.
void scatter_contact_block(const std::vector<std::size_t>& gdofs,
                           const std::vector<Real>& K, LinearSystem& sys) {
  const ConstraintTransform& tf = sys.transform;
  const int m = static_cast<int>(gdofs.size());
  for (int a = 0; a < m; ++a) {
    const DofExpansion& Ea = tf.expansion[gdofs[static_cast<std::size_t>(a)]];
    if (Ea.terms.empty() && Ea.constant == 0.0) continue;
    for (int b = 0; b < m; ++b) {
      const Real kab = K[static_cast<std::size_t>(a) * static_cast<std::size_t>(m) +
                         static_cast<std::size_t>(b)];
      if (kab == 0.0) continue;
      const DofExpansion& Eb = tf.expansion[gdofs[static_cast<std::size_t>(b)]];
      for (const DofTerm& ta : Ea.terms) {
        for (const DofTerm& tb : Eb.terms) {
          sys.rows.push_back(ta.eq);
          sys.cols.push_back(tb.eq);
          sys.vals.push_back(ta.coeff * kab * tb.coeff);
        }
        if (Eb.constant != 0.0)
          sys.rhs[static_cast<std::size_t>(ta.eq)] -= ta.coeff * kab * Eb.constant;
      }
    }
  }
}

// Broad phase (reference-mesh grid) + deformed-config projection: the master face nearest
// the deformed slave point `xs`, within the search distance. `have` is false when no
// candidate qualifies.
Projection nearest_master(const Mesh& mesh, const ResolvedContactPair& rp, const Vec3& xs,
                          const std::vector<Vec3>& u, bool& have) {
  Projection best;
  have = false;
  for (const int f : rp.grid.candidates(xs)) {
    const Projection pr =
        project_deformed(mesh, rp.master_faces[static_cast<std::size_t>(f)], xs, u);
    if (pr.dist > rp.grid.search_distance()) continue;
    if (!have || pr.dist < best.dist) {
      best = pr;
      have = true;
    }
  }
  return best;
}

// Consistent penalty tangent slope kt = -dF/dg at the current gap: kappa for the
// linear/hard law, the local exponential slope for the exponential law.
Real tangent_slope(const SurfaceBehavior& beh, Real kappa, Real gap) {
  if (beh.law == SurfaceBehavior::Law::Exponential && beh.c0 > 0.0)
    return kappa * std::exp((-gap) / beh.c0);
  return kappa;
}

// Add the normal contact force of one active slave node into the full-DOF internal force.
// f_int = kappa*g * dg/du = -Fmag * dg/du, with dg/du = +n on the slave and -N_i n on
// each master face node — the driver's residual r = f_ext - f_int then pushes the slave
// out and the reaction into the master.
void add_contact_force(std::size_t si, Real Fmag, const Vec3& n, const FaceEval& ev,
                       std::vector<Real>& f_int) {
  for (int d = 0; d < 3; ++d)
    f_int[si * kDofsPerNode + static_cast<std::size_t>(d)] -= Fmag * n[static_cast<std::size_t>(d)];
  for (std::size_t i = 0; i < ev.gnode.size(); ++i) {
    const std::size_t mi = static_cast<std::size_t>(ev.gnode[i]);
    for (int d = 0; d < 3; ++d)
      f_int[mi * kDofsPerNode + static_cast<std::size_t>(d)] += Fmag * ev.N[i] * n[static_cast<std::size_t>(d)];
  }
}

// Build the penalty tangent kt*(dg/du)⊗(dg/du) over the slave + master DOFs and scatter it
// through the constraint transform. `grad` = dg/du: +n on the slave, -N_i n on masters.
void add_contact_tangent(std::size_t si, Real kt, const Vec3& n, const FaceEval& ev,
                         LinearSystem& sys) {
  const int nf = static_cast<int>(ev.gnode.size());
  const int m = (1 + nf) * 3;
  std::vector<std::size_t> gdofs(static_cast<std::size_t>(m));
  std::vector<Real> grad(static_cast<std::size_t>(m), 0.0);
  for (int d = 0; d < 3; ++d) {
    gdofs[static_cast<std::size_t>(d)] = si * kDofsPerNode + static_cast<std::size_t>(d);
    grad[static_cast<std::size_t>(d)] = n[static_cast<std::size_t>(d)];
  }
  for (int i = 0; i < nf; ++i) {
    const std::size_t mi = static_cast<std::size_t>(ev.gnode[static_cast<std::size_t>(i)]);
    for (int d = 0; d < 3; ++d) {
      const std::size_t slot = static_cast<std::size_t>((1 + i) * 3 + d);
      gdofs[slot] = mi * kDofsPerNode + static_cast<std::size_t>(d);
      grad[slot] = -ev.N[static_cast<std::size_t>(i)] * n[static_cast<std::size_t>(d)];
    }
  }
  std::vector<Real> K(static_cast<std::size_t>(m) * static_cast<std::size_t>(m), 0.0);
  for (int a = 0; a < m; ++a)
    for (int b = 0; b < m; ++b)
      K[static_cast<std::size_t>(a) * static_cast<std::size_t>(m) + static_cast<std::size_t>(b)] =
          kt * grad[static_cast<std::size_t>(a)] * grad[static_cast<std::size_t>(b)];
  scatter_contact_block(gdofs, K, sys);
}

// Tangential relative displacement of the slave node vs the master footprint it projects
// onto, projected into the tangent plane (normal component removed). du_rel = u_slave -
// Σ N_i u_master_i; g_t = du_rel - (du_rel·n) n. Since the step starts from zero
// displacement, g_t IS the accumulated tangential slip within the increment (penalty
// stick reference), so the elastic predictor -stick*g_t needs no extra history.
Vec3 tangential_slip(std::size_t si, const Vec3& n, const FaceEval& ev,
                     const std::vector<Vec3>& u) {
  Vec3 um{0, 0, 0};
  for (std::size_t i = 0; i < ev.gnode.size(); ++i) {
    const std::size_t mi = static_cast<std::size_t>(ev.gnode[i]);
    for (int d = 0; d < 3; ++d)
      um[static_cast<std::size_t>(d)] += ev.N[i] * u[mi][static_cast<std::size_t>(d)];
  }
  Vec3 du{u[si][0] - um[0], u[si][1] - um[1], u[si][2] - um[2]};
  const Real dn = dot(du, n);
  return {du[0] - dn * n[0], du[1] - dn * n[1], du[2] - dn * n[2]};
}

// Add the tangential friction traction `t` of one active slave node into the full-DOF
// internal force. Mirrors add_contact_force with the normal replaced by the tangential
// traction: f_int gets -t on the slave and +N_i t on each master node (reaction), so the
// residual r = f_ext - f_int carries the friction force resisting the slave's tangential
// motion. `t` is the admissible (post-return-map) traction on the SLAVE.
void add_friction_force(std::size_t si, const Vec3& t, const FaceEval& ev,
                        std::vector<Real>& f_int) {
  for (int d = 0; d < 3; ++d)
    f_int[si * kDofsPerNode + static_cast<std::size_t>(d)] -= t[static_cast<std::size_t>(d)];
  for (std::size_t i = 0; i < ev.gnode.size(); ++i) {
    const std::size_t mi = static_cast<std::size_t>(ev.gnode[i]);
    for (int d = 0; d < 3; ++d)
      f_int[mi * kDofsPerNode + static_cast<std::size_t>(d)] += ev.N[i] * t[static_cast<std::size_t>(d)];
  }
}

// Add the consistent tangential friction tangent K = B Dmat Bᵀ across the slave + master
// DOFs, where B is the relative-displacement operator (row block +I on the slave, -N_i I
// on master face node i) and `Dmat` (3x3) is the tangent d(f_int,slave)/d(relative
// tangential displacement) of the friction traction. The two regimes supply different
// Dmat (both symmetric PSD-friendly, so the system stays direct-solve friendly):
//   STICK:  Dmat = stick * P            (elastic stick spring; P = I - n⊗n projector)
//   SLIP:   Dmat = (mu*p/|g_t|)(P - ŝ⊗ŝ)  (return-mapped: the force is capped at mu*p, so
//           the tangent is only the slip-DIRECTION change 1/|g_t|(P - ŝ⊗ŝ), NOT the stick
//           stiffness — this small tangent is what makes the slip Newton step converge).
// Dropping the normal-pressure cross term dp/du keeps Dmat symmetric (a standard symmetric
// penalty-friction linearization; the residual return-map still enforces the exact cone).
void add_friction_tangent(std::size_t si, const FaceEval& ev,
                          const std::array<std::array<Real, 3>, 3>& Dmat,
                          LinearSystem& sys) {
  const int nf = static_cast<int>(ev.gnode.size());
  const int m = (1 + nf) * 3;
  // Per-DOF B coefficient c_a (=+1 slave, -N_i master) and spatial direction dir(a).
  std::vector<std::size_t> gdofs(static_cast<std::size_t>(m));
  std::vector<Real> coeff(static_cast<std::size_t>(m), 0.0);
  std::vector<int> dir(static_cast<std::size_t>(m), 0);
  for (int d = 0; d < 3; ++d) {
    gdofs[static_cast<std::size_t>(d)] = si * kDofsPerNode + static_cast<std::size_t>(d);
    coeff[static_cast<std::size_t>(d)] = 1.0;
    dir[static_cast<std::size_t>(d)] = d;
  }
  for (int i = 0; i < nf; ++i) {
    const std::size_t mi = static_cast<std::size_t>(ev.gnode[static_cast<std::size_t>(i)]);
    for (int d = 0; d < 3; ++d) {
      const std::size_t slot = static_cast<std::size_t>((1 + i) * 3 + d);
      gdofs[slot] = mi * kDofsPerNode + static_cast<std::size_t>(d);
      coeff[slot] = -ev.N[static_cast<std::size_t>(i)];
      dir[slot] = d;
    }
  }
  std::vector<Real> K(static_cast<std::size_t>(m) * static_cast<std::size_t>(m), 0.0);
  for (int a = 0; a < m; ++a)
    for (int b = 0; b < m; ++b)
      K[static_cast<std::size_t>(a) * static_cast<std::size_t>(m) + static_cast<std::size_t>(b)] =
          coeff[static_cast<std::size_t>(a)] * coeff[static_cast<std::size_t>(b)] *
          Dmat[static_cast<std::size_t>(dir[static_cast<std::size_t>(a)])]
              [static_cast<std::size_t>(dir[static_cast<std::size_t>(b)])];
  scatter_contact_block(gdofs, K, sys);
}

// Tangent-plane projector P = I - n⊗n as a 3x3 matrix.
std::array<std::array<Real, 3>, 3> projector(const Vec3& n) {
  std::array<std::array<Real, 3>, 3> P{};
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      P[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] =
          (a == b ? 1.0 : 0.0) - n[static_cast<std::size_t>(a)] * n[static_cast<std::size_t>(b)];
  return P;
}

// The friction material tangent Dmat = d(f_int,slave)/d(g_t) for the current regime.
// STICK (not slipping): stick * P. SLIP: (mu*p/|g_t|)(P - ŝ⊗ŝ) with ŝ = g_t/|g_t| the slip
// direction — the small slip-direction tangent that keeps Newton convergent.
std::array<std::array<Real, 3>, 3> friction_tangent_matrix(bool slipping, Real stick,
                                                           Real mu, Real p_normal,
                                                           const Vec3& n, const Vec3& gt) {
  std::array<std::array<Real, 3>, 3> P = projector(n);
  if (!slipping) {  // stick: elastic stick spring stiffness on the tangent plane
    for (auto& row : P)
      for (Real& v : row) v *= stick;
    return P;
  }
  const Real gm = std::sqrt(dot(gt, gt));
  if (gm <= 0.0) return {};  // no slip direction yet -> no slip tangent
  const Real scale = mu * p_normal / gm;
  const Vec3 s{gt[0] / gm, gt[1] / gm, gt[2] / gm};  // unit slip direction
  std::array<std::array<Real, 3>, 3> D{};
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      D[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] =
          scale * (P[static_cast<std::size_t>(a)][static_cast<std::size_t>(b)] -
                   s[static_cast<std::size_t>(a)] * s[static_cast<std::size_t>(b)]);
  return D;
}

// Per-slave-node tributary contact area. For an element slave surface, sum the shape
// integral ∫N_i dA of each face over the nodes it touches (the consistent lumped area).
// For a node slave surface (no faces), fall back to the total master-face area split
// evenly over the slave nodes — the best area estimate without slave face topology. Keyed
// by dense node index; nodes not on a face get the fallback share.
std::unordered_map<Index, Real> slave_areas(const Mesh& mesh, const Surface& slave,
                                            const std::vector<Index>& slave_nodes,
                                            const std::vector<MasterFace>& mfaces) {
  std::unordered_map<Index, Real> area;
  if (slave.type == Surface::Type::Element) {
    for (const auto& [elem_id, face] : slave.faces) {
      const Index ei = mesh.element_index(elem_id);
      if (ei < 0) continue;
      const FaceSurface fs = face_surface_integrals(mesh, ei, face);
      for (std::size_t i = 0; i < fs.gnode.size(); ++i)
        area[fs.gnode[i]] += fs.load[i];  // ∫N_i dA (lumped tributary area)
    }
    return area;
  }
  // Node slave surface: split the total master area evenly over the slave nodes.
  Real total = 0.0;
  for (const MasterFace& mf : mfaces) {
    const FaceSurface fs = face_surface_integrals(mesh, mf.elem_index, mf.face);
    for (const Real l : fs.load) total += l;
  }
  const Real per = slave_nodes.empty() ? 0.0 : total / static_cast<Real>(slave_nodes.size());
  for (const Index n : slave_nodes) area[n] = per;
  return area;
}

}  // namespace

std::vector<ResolvedThermalContactPair> build_thermal_contact_pairs(const Model& model) {
  std::vector<ResolvedThermalContactPair> out;
  for (const ContactPair& pair : model.contact_pairs) {
    const GapConductance gc = model.contact_conductance(pair);
    const GapHeatGeneration gh = model.contact_heat_generation(pair);
    if (!gc.has && !gh.has) continue;  // no thermal-contact card on this pair
    if (pair.formulation != ContactPair::Formulation::NodeToSurface)
      throw std::runtime_error(
          "thermal contact is only implemented for TYPE=NODE TO SURFACE (mortar "
          "surface-to-surface is deferred; tasks.md 1.3)");
    const Surface* slave = model.mesh.surface(pair.slave_surface);
    const Surface* master = model.mesh.surface(pair.master_surface);
    if (slave == nullptr || master == nullptr)
      throw std::runtime_error("*CONTACT PAIR references unknown surface for thermal contact");
    std::vector<MasterFace> faces = expand_master_faces(model.mesh, *master);
    if (faces.empty())
      throw std::runtime_error("thermal *CONTACT PAIR master surface has no element faces");
    std::vector<Index> snodes = resolve_slave_nodes(model.mesh, *slave);
    Real search = pair.search;
    if (search <= 0.0)
      search = 2.0 * element_size(
                         model.mesh,
                         model.mesh.elements()[static_cast<std::size_t>(faces.front().elem_index)]);
    const std::unordered_map<Index, Real> amap =
        slave_areas(model.mesh, *slave, snodes, faces);
    std::vector<Real> areas;
    areas.reserve(snodes.size());
    for (const Index n : snodes) {
      const auto it = amap.find(n);
      areas.push_back(it != amap.end() ? it->second : 0.0);
    }
    const Real adjust = 1e-3 * element_size(
                                   model.mesh,
                                   model.mesh.elements()[static_cast<std::size_t>(
                                       faces.front().elem_index)]);
    ResolvedThermalContactPair tp{snodes,
                                  std::move(areas),
                                  std::move(faces),
                                  MasterGrid(model.mesh, {}, search),
                                  gc.h,
                                  gh.q,
                                  adjust};
    // Rebuild the grid over the struct's OWN master_faces (the grid stores a pointer to
    // the face vector, so it must reference the persistent member, not the local temp).
    tp.grid = MasterGrid(model.mesh, tp.master_faces, search);
    out.push_back(std::move(tp));
  }
  return out;
}

int add_thermal_contact(const Model& model,
                        const std::vector<ResolvedThermalContactPair>& pairs,
                        const std::vector<Vec3>& u, std::vector<Index>& k_rows,
                        std::vector<Index>& k_cols, std::vector<Real>& k_vals,
                        std::vector<Real>& source) {
  const Mesh& mesh = model.mesh;
  int active = 0;
  for (const ResolvedThermalContactPair& tp : pairs) {
    for (std::size_t s = 0; s < tp.slave_nodes.size(); ++s) {
      const std::size_t si = static_cast<std::size_t>(tp.slave_nodes[s]);
      const Vec3 xs{mesh.nodes()[si].x[0] + u[si][0], mesh.nodes()[si].x[1] + u[si][1],
                    mesh.nodes()[si].x[2] + u[si][2]};
      // Nearest master face in the current (deformed) config; skip if none in contact.
      bool have = false;
      Projection best;
      for (const int fi : tp.grid.candidates(xs)) {
        const Projection pr =
            project_deformed(mesh, tp.master_faces[static_cast<std::size_t>(fi)], xs, u);
        if (pr.dist > tp.grid.search_distance()) continue;
        if (!have || pr.dist < best.dist) { best = pr; have = true; }
      }
      if (!have || best.gap >= tp.adjust) continue;  // not in contact -> no gap heat flow
      ++active;
      const Real A = tp.slave_area[s];
      if (A <= 0.0) continue;

      // Conductance block: Q = h A (T_s - Σ N_i T_mi). Add h A on (s,s), -h A N_i on
      // (s, mi) and (mi, s), and +h A N_i N_j on (mi, mj) — the symmetric conductance
      // operator of the scalar coupling B = [1, -N_i], K = h A B^T B.
      const Real hA = tp.h * A;
      const FaceEval& ev = best.ev;
      const int nf = static_cast<int>(ev.gnode.size());
      if (hA > 0.0) {
        k_rows.push_back(static_cast<Index>(si)); k_cols.push_back(static_cast<Index>(si));
        k_vals.push_back(hA);
        for (int i = 0; i < nf; ++i) {
          const Index mi = ev.gnode[static_cast<std::size_t>(i)];
          const Real Ni = ev.N[static_cast<std::size_t>(i)];
          k_rows.push_back(static_cast<Index>(si)); k_cols.push_back(mi);
          k_vals.push_back(-hA * Ni);
          k_rows.push_back(mi); k_cols.push_back(static_cast<Index>(si));
          k_vals.push_back(-hA * Ni);
          for (int j = 0; j < nf; ++j) {
            const Index mj = ev.gnode[static_cast<std::size_t>(j)];
            const Real Nj = ev.N[static_cast<std::size_t>(j)];
            k_rows.push_back(mi); k_cols.push_back(mj);
            k_vals.push_back(hA * Ni * Nj);
          }
        }
      }
      // Gap heat generation: total q_gap*A generated, split evenly between the surfaces
      // — half to the slave node, half spread over the master face by N_i.
      if (tp.q_gap != 0.0) {
        const Real half = 0.5 * tp.q_gap * A;
        source[si] += half;
        for (int i = 0; i < nf; ++i)
          source[static_cast<std::size_t>(ev.gnode[static_cast<std::size_t>(i)])] +=
              half * ev.N[static_cast<std::size_t>(i)];
      }
    }
  }
  return active;
}

int add_contact(const Model& model, const std::vector<ResolvedContactPair>& pairs,
                const std::vector<Vec3>& u, LinearSystem& sys,
                std::vector<Real>& f_int) {
  const Mesh& mesh = model.mesh;
  int active = 0;
  for (const ResolvedContactPair& rp : pairs) {
    for (const Index sni : rp.slave_nodes) {
      const std::size_t si = static_cast<std::size_t>(sni);
      const Vec3 xs{mesh.nodes()[si].x[0] + u[si][0], mesh.nodes()[si].x[1] + u[si][1],
                    mesh.nodes()[si].x[2] + u[si][2]};
      bool have = false;
      Projection best = nearest_master(mesh, rp, xs, u, have);
      // Apply the *CLEARANCE offset: the operator sees g_eff = g_geometric + clearance,
      // so a positive clearance opens the interface regardless of the as-meshed geometry.
      best.gap += rp.clearance;
      // Active set: penetrating (gap < 0), or within the near-contact adjust band
      // (0 <= gap < adjust) where the penalty stiffness kills the rigid mode but the force
      // is still zero. A clearance beyond the band releases the pair (tensile) -> inactive.
      if (!have || best.gap >= rp.adjust) continue;
      if (best.gap < 0.0) ++active;  // count genuinely-penetrating pairs

      const Real Fmag = contact_force_weight(rp.behavior, rp.penalty, best.gap);
      const Real kt = tangent_slope(rp.behavior, rp.penalty, best.gap);
      add_contact_force(si, Fmag, best.normal, best.ev, f_int);
      add_contact_tangent(si, kt, best.normal, best.ev, sys);

      // Tangential Coulomb friction (only under genuine normal pressure, gap < 0). The
      // elastic stick predictor t_trial = -stick * g_t is return-mapped onto the friction
      // cone |t| <= mu*Fmag: inside -> stick (t = t_trial), outside -> slip (t capped at
      // mu*Fmag along g_t). Both regimes add the traction to the residual and the stick
      // tangent to the system. (spec: contact — tangential contact / stick-slip.)
      if (rp.friction.has && rp.friction.mu > 0.0 && best.gap < 0.0) {
        const Vec3 gt = tangential_slip(si, best.normal, best.ev, u);
        const Vec3 t_trial{-rp.stick * gt[0], -rp.stick * gt[1], -rp.stick * gt[2]};
        bool slipping = false;
        const Vec3 t = friction_return_map(t_trial, rp.friction.mu, Fmag, slipping);
        add_friction_force(si, t, best.ev, f_int);
        // Consistent regime tangent: stick spring (stick*P) while stuck, the small
        // slip-direction tangent once slipping — the slip tangent is what lets the Newton
        // step converge instead of chattering against the capped force.
        const std::array<std::array<Real, 3>, 3> Dmat = friction_tangent_matrix(
            slipping, rp.stick, rp.friction.mu, Fmag, best.normal, gt);
        add_friction_tangent(si, best.ev, Dmat, sys);
      }
    }
  }
  return active;
}

std::vector<ContactPoint> recover_contact(const Model& model,
                                          const std::vector<ResolvedContactPair>& pairs,
                                          const std::vector<Vec3>& u) {
  const Mesh& mesh = model.mesh;
  std::vector<ContactPoint> out;
  for (const ResolvedContactPair& rp : pairs) {
    // Tributary area per slave node ~ total master-face area / number of slave nodes, so
    // the reported pressure p = F/area is the physical interface pressure.
    Real master_area = 0.0;
    for (const MasterFace& mf : rp.master_faces) {
      const FaceSurface fs = face_surface_integrals(mesh, mf.elem_index, mf.face);
      for (const Real l : fs.load) master_area += l;
    }
    const Real area = rp.slave_nodes.empty()
                          ? 0.0
                          : master_area / static_cast<Real>(rp.slave_nodes.size());

    for (const Index sni : rp.slave_nodes) {
      const std::size_t si = static_cast<std::size_t>(sni);
      const Vec3 xs{mesh.nodes()[si].x[0] + u[si][0], mesh.nodes()[si].x[1] + u[si][1],
                    mesh.nodes()[si].x[2] + u[si][2]};
      bool have = false;
      Projection best = nearest_master(mesh, rp, xs, u, have);
      best.gap += rp.clearance;  // *CLEARANCE offset (same as add_contact)

      ContactPoint cp;
      cp.node_id = mesh.nodes()[si].id;
      // OPEN when no master found or the node is clear beyond the adjust band; CLOSED when
      // penetrating (the active set the operator uses to build a contact force).
      cp.closed = have && best.gap < 0.0;
      cp.gap = have ? best.gap : 0.0;
      if (cp.closed) {
        const Real Fmag = contact_force_weight(rp.behavior, rp.penalty, best.gap);
        cp.p = area > 0.0 ? Fmag / area : 0.0;  // normal contact pressure
        if (rp.friction.has && rp.friction.mu > 0.0) {
          const Vec3 gt = tangential_slip(si, best.normal, best.ev, u);
          const Vec3 t_trial{-rp.stick * gt[0], -rp.stick * gt[1], -rp.stick * gt[2]};
          bool slipping = false;
          const Vec3 t = friction_return_map(t_trial, rp.friction.mu, Fmag, slipping);
          const Real tmag = std::sqrt(dot(t, t));
          cp.tau = area > 0.0 ? tmag / area : 0.0;  // tangential traction (per area)
        }
      }
      out.push_back(cp);
    }
  }
  return out;
}

}  // namespace cxpp::fem
