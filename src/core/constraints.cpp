#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

#include "calculixpp/core/model.hpp"

// Expansion of the high-level constraint cards (*MPC / *RIGID BODY / *COUPLING /
// *TIE) into flat linear Equations (spec: constraints). Every constraint type
// reduces to  sum_i coeff_i u_(node_i, comp_i) = 0  with the first term dependent;
// the assembly then eliminates each dependent DOF. Small-displacement kinematics:
// rigid links and couplings tie translations directly (no finite rotation), which is
// the linear-static regime these tasks target.

namespace cxpp {
namespace {

// Node coordinate lookup (throws on an unknown node).
Vec3 node_x(const Mesh& mesh, Index node_id) {
  const Index ni = mesh.node_index(node_id);
  if (ni < 0)
    throw std::runtime_error("constraint references unknown node " +
                             std::to_string(node_id));
  return mesh.nodes()[static_cast<std::size_t>(ni)].x;
}

// Resolve a node-set name (or, if empty, an explicit id list) to node ids.
std::vector<Index> resolve_nodes(const Mesh& mesh, const std::string& nset,
                                 const std::vector<Index>& explicit_ids) {
  if (!nset.empty()) {
    const std::vector<Index>* ids = mesh.nset(nset);
    if (ids == nullptr)
      throw std::runtime_error("constraint references unknown node set '" + nset + "'");
    return *ids;
  }
  return explicit_ids;
}

// Tie a whole translation (all 3 comps) of `slave` to `master`:  u_slave - u_master = 0
// per component (slave dependent). Used by BEAM MPCs, rigid bodies (no rotation),
// kinematic couplings, and matching ties.
void tie_translation(std::vector<Equation>& out, Index slave, Index master,
                     const char* origin) {
  for (int c = 1; c <= 3; ++c)
    out.push_back(Equation{{{slave, c, 1.0}, {master, c, -1.0}}, origin});
}

// STRAIGHT MPC: keep node p on the straight line through a (first) and b (second).
// Small-displacement linearization: the transverse displacement of p relative to the
// a-b line is zero. Build two equations from the two unit vectors t1,t2 spanning the
// plane normal to the line direction e:  t_k . (u_p - u_a) - s (t_k.(u_b - u_a)) = 0,
// where s = |x_p - x_a| / |x_b - x_a| is the along-line parameter of p. p's leading
// (dependent) DOF is the component of t_k with the largest magnitude.
void straight_mpc(std::vector<Equation>& out, const Mesh& mesh,
                  const std::vector<Index>& nodes) {
  if (nodes.size() < 3) return;  // need a, b, and at least one governed node
  const Vec3 xa = node_x(mesh, nodes[0]);
  const Vec3 xb = node_x(mesh, nodes[1]);
  std::array<Real, 3> e{xb[0] - xa[0], xb[1] - xa[1], xb[2] - xa[2]};
  const Real len = std::sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
  if (len <= 0.0) throw std::runtime_error("*MPC STRAIGHT has coincident base nodes");
  for (Real& v : e) v /= len;
  // Two orthonormal vectors t1,t2 perpendicular to e.
  const std::array<Real, 3> ref =
      std::abs(e[0]) < 0.9 ? std::array<Real, 3>{1, 0, 0} : std::array<Real, 3>{0, 1, 0};
  std::array<Real, 3> t1{e[1] * ref[2] - e[2] * ref[1], e[2] * ref[0] - e[0] * ref[2],
                         e[0] * ref[1] - e[1] * ref[0]};
  const Real n1 = std::sqrt(t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]);
  for (Real& v : t1) v /= n1;
  std::array<Real, 3> t2{e[1] * t1[2] - e[2] * t1[1], e[2] * t1[0] - e[0] * t1[2],
                         e[0] * t1[1] - e[1] * t1[0]};
  const auto emit = [&](Index p, const std::array<Real, 3>& t) {
    const Vec3 xp = node_x(mesh, p);
    std::array<Real, 3> d{xp[0] - xa[0], xp[1] - xa[1], xp[2] - xa[2]};
    const Real s = (d[0] * e[0] + d[1] * e[1] + d[2] * e[2]) / len;
    // t.(u_p) - (1-s) t.(u_a) - s t.(u_b) = 0. Dependent term: p's largest |t| comp.
    int dep = 0;
    for (int c = 1; c < 3; ++c)
      if (std::abs(t[static_cast<std::size_t>(c)]) > std::abs(t[static_cast<std::size_t>(dep)])) dep = c;
    Equation eq;
    eq.origin = "*MPC STRAIGHT";
    eq.terms.push_back({p, dep + 1, t[static_cast<std::size_t>(dep)]});
    for (int c = 0; c < 3; ++c)
      if (c != dep && t[static_cast<std::size_t>(c)] != 0.0)
        eq.terms.push_back({p, c + 1, t[static_cast<std::size_t>(c)]});
    for (int c = 0; c < 3; ++c) {
      if (t[static_cast<std::size_t>(c)] == 0.0) continue;
      eq.terms.push_back({nodes[0], c + 1, -(1.0 - s) * t[static_cast<std::size_t>(c)]});
      eq.terms.push_back({nodes[1], c + 1, -s * t[static_cast<std::size_t>(c)]});
    }
    out.push_back(std::move(eq));
  };
  for (std::size_t i = 2; i < nodes.size(); ++i) {
    emit(nodes[i], t1);
    emit(nodes[i], t2);
  }
}

// PLANE MPC: keep each governed node in the plane through the first three nodes.
// Small-displacement: the out-of-plane displacement of node p (relative to the plane
// spanned by nodes a,b,c) is zero. Normal n = (xb-xa) x (xc-xa). Equation:
// n.(u_p) - w_a n.u_a - w_b n.u_b - w_c n.u_c = 0 with barycentric weights of p in the
// plane a,b,c. For simplicity the base-node coupling uses the mean (w=1/3), which is
// exact when p lies at the centroid and is a consistent rigid-plane constraint
// otherwise (the plane translates/rotates with a,b,c). Dependent DOF: largest |n|.
void plane_mpc(std::vector<Equation>& out, const Mesh& mesh,
               const std::vector<Index>& nodes) {
  if (nodes.size() < 4) return;  // 3 base nodes define the plane, need >=1 governed
  const Vec3 xa = node_x(mesh, nodes[0]);
  const Vec3 xb = node_x(mesh, nodes[1]);
  const Vec3 xc = node_x(mesh, nodes[2]);
  std::array<Real, 3> u{xb[0] - xa[0], xb[1] - xa[1], xb[2] - xa[2]};
  std::array<Real, 3> v{xc[0] - xa[0], xc[1] - xa[1], xc[2] - xa[2]};
  std::array<Real, 3> n{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2],
                        u[0] * v[1] - u[1] * v[0]};
  const Real nn = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  if (nn <= 0.0) throw std::runtime_error("*MPC PLANE base nodes are collinear");
  for (Real& x : n) x /= nn;
  int dep = 0;
  for (int c = 1; c < 3; ++c)
    if (std::abs(n[static_cast<std::size_t>(c)]) > std::abs(n[static_cast<std::size_t>(dep)])) dep = c;
  for (std::size_t i = 3; i < nodes.size(); ++i) {
    Equation eq;
    eq.origin = "*MPC PLANE";
    eq.terms.push_back({nodes[i], dep + 1, n[static_cast<std::size_t>(dep)]});
    for (int c = 0; c < 3; ++c)
      if (c != dep && n[static_cast<std::size_t>(c)] != 0.0)
        eq.terms.push_back({nodes[i], c + 1, n[static_cast<std::size_t>(c)]});
    for (int b = 0; b < 3; ++b)
      for (int c = 0; c < 3; ++c)
        if (n[static_cast<std::size_t>(c)] != 0.0)
          eq.terms.push_back({nodes[static_cast<std::size_t>(b)], c + 1,
                              -(1.0 / 3.0) * n[static_cast<std::size_t>(c)]});
    out.push_back(std::move(eq));
  }
}

// *RIGID BODY: each set node translation = reference translation + omega x (x - x_ref).
// Without a rotation node the rotational part is dropped (pure translation tie). With
// a rotation node, the three rotation DOFs live on that node's comps 1..3 and the
// cross product is linearized (small rotation). Equation per set node per comp:
//   u_node,i - u_ref,i - sum_j eps_{ijk} r_k theta_j? ...  small-rotation form:
//   u_node = u_ref + theta x r,  r = x_node - x_ref.
void rigid_body(std::vector<Equation>& out, const Mesh& mesh, const RigidBody& rb) {
  const std::vector<Index> nodes = resolve_nodes(mesh, rb.nset, rb.nodes);
  const Vec3 xref = node_x(mesh, rb.ref_node);
  for (const Index nd : nodes) {
    if (nd == rb.ref_node) continue;
    if (rb.rot_node == 0) {
      tie_translation(out, nd, rb.ref_node, "*RIGID BODY");
      continue;
    }
    const Vec3 xn = node_x(mesh, nd);
    const std::array<Real, 3> r{xn[0] - xref[0], xn[1] - xref[1], xn[2] - xref[2]};
    // u_node,i = u_ref,i + (theta x r)_i. (theta x r) = (t2 r3 - t3 r2, t3 r1 - t1 r3,
    // t1 r2 - t2 r1). theta components live on rot_node comps 1..3.
    const std::array<std::array<Real, 3>, 3> cross{{{0, -r[2], r[1]},
                                                    {r[2], 0, -r[0]},
                                                    {-r[1], r[0], 0}}};
    for (int i = 0; i < 3; ++i) {
      Equation eq;
      eq.origin = "*RIGID BODY";
      eq.terms.push_back({nd, i + 1, 1.0});
      eq.terms.push_back({rb.ref_node, i + 1, -1.0});
      for (int j = 0; j < 3; ++j)
        if (cross[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] != 0.0)
          eq.terms.push_back({rb.rot_node, j + 1,
                              -cross[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]});
      out.push_back(std::move(eq));
    }
  }
}

// *COUPLING. KINEMATIC: every surface node's selected DOFs equal the reference node's
// (rigid translation tie). DISTRIBUTING: the reference node motion is the weighted
// average of the surface nodes -> one averaging equation per DOF with the reference
// node dependent (equal weights). (spec: constraints — couplings.)
void coupling(std::vector<Equation>& out, const Coupling& cp) {
  if (cp.kind == Coupling::Kind::Kinematic) {
    for (const Index nd : cp.nodes) {
      if (nd == cp.ref_node) continue;
      for (const int c : cp.dofs)
        out.push_back(Equation{{{nd, c, 1.0}, {cp.ref_node, c, -1.0}}, "*COUPLING KINEMATIC"});
    }
    return;
  }
  // Distributing: u_ref,c = mean(u_surface,c). Reference DOF is dependent.
  const Real w = cp.nodes.empty() ? 0.0 : 1.0 / static_cast<Real>(cp.nodes.size());
  for (const int c : cp.dofs) {
    Equation eq;
    eq.origin = "*DISTRIBUTING COUPLING";
    eq.terms.push_back({cp.ref_node, c, 1.0});
    for (const Index nd : cp.nodes) eq.terms.push_back({nd, c, -w});
    out.push_back(std::move(eq));
  }
}

// *TIE for matching surfaces: pair each slave node with the coincident master node
// (nearest within tolerance) and tie all three translations. Non-matching (mortar)
// surfaces are not fully supported: a slave node with no master within tolerance is
// skipped (partial — see tasks.md 5.5).
void tie(std::vector<Equation>& out, const Mesh& mesh, const Tie& t) {
  const Real tol = t.tolerance > 0.0 ? t.tolerance : 1e-6;
  for (const Index sn : t.slave_nodes) {
    const Vec3 xs = node_x(mesh, sn);
    Index best = -1;
    Real best_d2 = tol * tol;
    for (const Index mn : t.master_nodes) {
      if (mn == sn) continue;
      const Vec3 xm = node_x(mesh, mn);
      const Real d2 = (xs[0] - xm[0]) * (xs[0] - xm[0]) +
                      (xs[1] - xm[1]) * (xs[1] - xm[1]) +
                      (xs[2] - xm[2]) * (xs[2] - xm[2]);
      if (d2 <= best_d2) { best_d2 = d2; best = mn; }
    }
    if (best >= 0) tie_translation(out, sn, best, "*TIE");
  }
}

}  // namespace

std::vector<Equation> Model::expand_constraints() const {
  std::vector<Equation> out = equations;  // raw *EQUATION relations first
  for (const Mpc& m : mpcs) {
    if (m.kind == Mpc::Kind::Beam) {
      if (m.nodes.size() >= 2) tie_translation(out, m.nodes[0], m.nodes[1], "*MPC BEAM");
    } else if (m.kind == Mpc::Kind::Straight) {
      straight_mpc(out, mesh, m.nodes);
    } else if (m.kind == Mpc::Kind::Plane) {
      plane_mpc(out, mesh, m.nodes);
    }
    // Kind::User: hook not registered -> no equations (documented, see tasks 5.2).
  }
  for (const RigidBody& rb : rigid_bodies) rigid_body(out, mesh, rb);
  for (const Coupling& cp : couplings) coupling(out, cp);
  for (const Tie& t : ties) tie(out, mesh, t);
  return out;
}

}  // namespace cxpp
