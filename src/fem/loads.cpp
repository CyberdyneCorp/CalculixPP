#include "calculixpp/fem/loads.hpp"

#include <array>
#include <cmath>
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

// 3-point order-2 triangle Gauss rule (reference triangle area 1/2).
struct TriGP {
  Real xi, eta, w;
};
constexpr std::array<TriGP, 3> kTri{{
    {1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0},
    {2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0},
    {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0},
}};

// Face shape functions and derivatives (T3 for 3 nodes, T6 for 6 nodes).
void face_shape(int nf, Real xi, Real eta, std::array<Real, 6>& N,
                std::array<Real, 6>& dxi, std::array<Real, 6>& deta) {
  const Real l1 = 1.0 - xi - eta, l2 = xi, l3 = eta;
  if (nf == 3) {
    N = {l1, l2, l3, 0, 0, 0};
    dxi = {-1, 1, 0, 0, 0, 0};
    deta = {-1, 0, 1, 0, 0, 0};
    return;
  }
  N = {l1 * (2 * l1 - 1), l2 * (2 * l2 - 1), l3 * (2 * l3 - 1),
       4 * l1 * l2,       4 * l2 * l3,       4 * l3 * l1};
  dxi = {-(4 * l1 - 1), 4 * l2 - 1, 0, 4 * (l1 - l2), 4 * l3, -4 * l3};
  deta = {-(4 * l1 - 1), 0, 4 * l3 - 1, -4 * l2, 4 * l2, 4 * (l1 - l3)};
}

// Accumulate consistent nodal forces for a pressure `p` on element face (*DLOAD
// P<face> or *DSLOAD), F_i = -p * integral( N_i * (t_xi x t_eta) ), into f.
void accumulate_pressure_face(const Mesh& mesh, Index ei, int face, Real p,
                              std::vector<Real>& f) {
  const Element& elem = mesh.elements()[static_cast<std::size_t>(ei)];
  const std::vector<int> fn = face_nodes(elem.type, face);
  const int nf = static_cast<int>(fn.size());

  std::vector<Vec3> x(static_cast<std::size_t>(nf));
  std::vector<Index> gnode(static_cast<std::size_t>(nf));
  for (int i = 0; i < nf; ++i) {
    const Index ni = mesh.node_index(elem.nodes[static_cast<std::size_t>(fn[static_cast<std::size_t>(i)])]);
    gnode[static_cast<std::size_t>(i)] = ni;
    x[static_cast<std::size_t>(i)] = mesh.nodes()[static_cast<std::size_t>(ni)].x;
  }

  std::array<Real, 6> N{}, dxi{}, deta{};
  std::vector<Vec3> force(static_cast<std::size_t>(nf), Vec3{0, 0, 0});
  for (const TriGP& gp : kTri) {
    face_shape(nf, gp.xi, gp.eta, N, dxi, deta);
    Vec3 txi{0, 0, 0}, teta{0, 0, 0};
    for (int i = 0; i < nf; ++i)
      for (int d = 0; d < 3; ++d) {
        txi[static_cast<std::size_t>(d)] += dxi[static_cast<std::size_t>(i)] * x[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)];
        teta[static_cast<std::size_t>(d)] += deta[static_cast<std::size_t>(i)] * x[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)];
      }
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
      f[static_cast<std::size_t>(gnode[static_cast<std::size_t>(i)]) * kDofsPerNode + static_cast<std::size_t>(d)] +=
          force[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)];
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

  std::array<std::array<Real, 3>, kMaxNodes> grad{};
  const auto do_element = [&](std::size_t e) {
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

}  // namespace

std::vector<int> face_nodes(ElementType type, int face) {
  // Pressure faces are currently implemented for tetrahedra only. Hex/wedge face
  // topology (6 quad + optional tri faces, different node counts) is not wired for
  // *DLOAD P<face>/*DSLOAD yet; body loads (GRAV/CENTRIF) work for all topologies.
  if (type != ElementType::C3D4 && type != ElementType::C3D10)
    throw std::runtime_error(
        "pressure faces (*DLOAD P<face>) are supported only for C3D4/C3D10 "
        "tetrahedra; hex/wedge pressure faces are not implemented yet");
  if (face < 1 || face > 4) throw std::runtime_error("tet face must be 1..4");
  const int nf = (type == ElementType::C3D4) ? 3 : 6;
  std::vector<int> out(static_cast<std::size_t>(nf));
  for (int i = 0; i < nf; ++i) out[static_cast<std::size_t>(i)] = kFace[face - 1][i];
  return out;
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

  for (const Dload& dl : model.dloads) {
    const Index ei = mesh.element_index(dl.elem_id);
    if (ei < 0) throw std::runtime_error("*DLOAD references unknown element");
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
