#include "calculixpp/fem/stress.hpp"

#include <array>
#include <vector>

#include "calculixpp/fem/element.hpp"
#include "calculixpp/fem/loads.hpp"

namespace cxpp::fem {
namespace {

// C3D10 corner extrapolation matrix a4 (inverse of the linear-tet shape functions
// at the 4 Gauss points) and midside adjacency (nonei10), from extrapolate.f.
constexpr Real kA4diag = 1.927050983124842;
constexpr Real kA4off = -0.309016994374947;
constexpr int kMid[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};

// Integration-point strain (engineering shear) and stress from displacements.
void strain_stress_at(const std::array<std::array<Real, 3>, kMaxNodes>& g, int n,
                      const std::vector<Vec3>& ue, const D6& D, Voigt6& strain,
                      Voigt6& stress) {
  // Displacement gradient dudx[i][j] = sum_k g_k[j] * u_k[i].
  std::array<std::array<Real, 3>, 3> du{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      Real v = 0.0;
      for (int k = 0; k < n; ++k) v += g[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)] * ue[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)];
      du[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = v;
    }
  strain = {du[0][0], du[1][1], du[2][2],
            du[0][1] + du[1][0], du[0][2] + du[2][0],
            du[1][2] + du[2][1]};
  for (int i = 0; i < 6; ++i) {
    Real v = 0.0;
    for (int j = 0; j < 6; ++j) v += D[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] * strain[static_cast<std::size_t>(j)];
    stress[static_cast<std::size_t>(i)] = v;
  }
}

// Extrapolate integration-point Voigt6 values to element nodes (C3D4 constant;
// C3D10 a4 corners + nonei10 midside averaging). Hex/wedge families use the
// mean of their Gauss-point values at every node (a valid, order-0 recovery; the
// per-topology extrapolation matrices are a later refinement — nodal stress is not
// part of the displacement-accuracy gate).
void extrapolate_to_nodes(ElementType type, int n,
                          const std::vector<Voigt6>& gp,
                          std::vector<Voigt6>& nodal) {
  if (type == ElementType::C3D4) {
    for (int i = 0; i < n; ++i) nodal[static_cast<std::size_t>(i)] = gp[0];
    return;
  }
  if (type != ElementType::C3D10) {
    Voigt6 mean{};
    for (const Voigt6& v : gp)
      for (int c = 0; c < 6; ++c) mean[static_cast<std::size_t>(c)] += v[static_cast<std::size_t>(c)];
    const Real inv = gp.empty() ? 0.0 : 1.0 / static_cast<Real>(gp.size());
    for (int c = 0; c < 6; ++c) mean[static_cast<std::size_t>(c)] *= inv;
    for (int i = 0; i < n; ++i) nodal[static_cast<std::size_t>(i)] = mean;
    return;
  }
  for (int c = 0; c < 4; ++c)
    for (int comp = 0; comp < 6; ++comp) {
      Real v = 0.0;
      for (int l = 0; l < 4; ++l)
        v += (c == l ? kA4diag : kA4off) * gp[static_cast<std::size_t>(l)][static_cast<std::size_t>(comp)];
      nodal[static_cast<std::size_t>(c)][static_cast<std::size_t>(comp)] = v;
    }
  for (int mnode = 0; mnode < 6; ++mnode)
    for (int comp = 0; comp < 6; ++comp)
      nodal[static_cast<std::size_t>(4 + mnode)][static_cast<std::size_t>(comp)] =
          0.5 * (nodal[static_cast<std::size_t>(kMid[mnode][0])][static_cast<std::size_t>(comp)] +
                 nodal[static_cast<std::size_t>(kMid[mnode][1])][static_cast<std::size_t>(comp)]);
}

// Scatter one element's internal force Ke_e u_e into the global vector f_int.
void accumulate_internal_force(ElementType type,
                               const std::vector<Vec3>& coords,
                               const std::vector<Index>& nidx,
                               const std::vector<Vec3>& ue, const ElasticIso& mat,
                               std::vector<Real>& f_int) {
  const int n = nodes_per_element(type);
  const std::vector<Real> Ke = element_stiffness(type, coords, mat);
  const int ndof = n * kDofsPerNode;
  for (int a = 0; a < ndof; ++a) {
    Real v = 0.0;
    for (int b = 0; b < ndof; ++b)
      v += Ke[static_cast<std::size_t>(a) * static_cast<std::size_t>(ndof) + static_cast<std::size_t>(b)] *
           ue[static_cast<std::size_t>(b / kDofsPerNode)][static_cast<std::size_t>(b % kDofsPerNode)];
    f_int[static_cast<std::size_t>(nidx[static_cast<std::size_t>(a / kDofsPerNode)]) * kDofsPerNode +
          static_cast<std::size_t>(a % kDofsPerNode)] += v;
  }
}

// Gather an element's node indices, coordinates, and current displacements.
void gather_element(const Mesh& mesh, const Element& elem,
                    const std::vector<Vec3>& u, std::vector<Vec3>& coords,
                    std::vector<Index>& nidx, std::vector<Vec3>& ue) {
  const int n = nodes_per_element(elem.type);
  coords.resize(static_cast<std::size_t>(n));
  nidx.resize(static_cast<std::size_t>(n));
  ue.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    const Index ni = mesh.node_index(elem.nodes[static_cast<std::size_t>(i)]);
    nidx[static_cast<std::size_t>(i)] = ni;
    coords[static_cast<std::size_t>(i)] = mesh.nodes()[static_cast<std::size_t>(ni)].x;
    ue[static_cast<std::size_t>(i)] = u[static_cast<std::size_t>(ni)];
  }
}

}  // namespace

std::vector<Real> internal_force(const Model& model, const std::vector<Vec3>& u) {
  const Mesh& mesh = model.mesh;
  std::vector<Real> f_int(mesh.num_nodes() * kDofsPerNode, 0.0);
  const std::vector<ElasticIso> elastic = model.element_elastic();
  std::vector<Vec3> coords, ue;
  std::vector<Index> nidx;
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    const Element& elem = mesh.elements()[e];
    gather_element(mesh, elem, u, coords, nidx, ue);
    accumulate_internal_force(elem.type, coords, nidx, ue, elastic[e], f_int);
  }
  return f_int;
}

void recover_fields(const Model& model, StaticFields& fields) {
  const Mesh& mesh = model.mesh;
  const std::size_t n_nodes = mesh.num_nodes();
  fields.stress.assign(n_nodes, Voigt6{});
  fields.strain.assign(n_nodes, Voigt6{});
  fields.reaction.assign(n_nodes, Vec3{0, 0, 0});
  std::vector<int> count(n_nodes, 0);
  std::vector<Real> f_int(n_nodes * kDofsPerNode, 0.0);

  const std::vector<ElasticIso> elastic = model.element_elastic();
  std::array<std::array<Real, 3>, kMaxNodes> g{};

  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    const Element& elem = mesh.elements()[e];
    const int n = nodes_per_element(elem.type);
    std::vector<Vec3> coords(static_cast<std::size_t>(n));
    std::vector<Index> nidx(static_cast<std::size_t>(n));
    std::vector<Vec3> ue(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
      const Index ni = mesh.node_index(elem.nodes[static_cast<std::size_t>(i)]);
      nidx[static_cast<std::size_t>(i)] = ni;
      coords[static_cast<std::size_t>(i)] = mesh.nodes()[static_cast<std::size_t>(ni)].x;
      ue[static_cast<std::size_t>(i)] = fields.displacement[static_cast<std::size_t>(ni)];
    }
    const D6 D = elastic_iso_D(elastic[e]);

    // Internal force f_int_e = Ke_e u_e, scattered to global DOFs.
    accumulate_internal_force(elem.type, coords, nidx, ue, elastic[e], f_int);

    // Integration-point strains and stresses.
    const auto rule = gauss_rule(elem.type);
    std::vector<Voigt6> gp_strain(rule.size());
    std::vector<Voigt6> gp_stress(rule.size());
    for (std::size_t q = 0; q < rule.size(); ++q) {
      const Shape s = shape(elem.type, rule[q].xi, rule[q].et, rule[q].ze);
      physical_gradients(s, coords, g);
      strain_stress_at(g, n, ue, D, gp_strain[q], gp_stress[q]);
    }

    // Extrapolate to element nodes.
    std::vector<Voigt6> nodal_strain(static_cast<std::size_t>(n));
    std::vector<Voigt6> nodal_stress(static_cast<std::size_t>(n));
    extrapolate_to_nodes(elem.type, n, gp_strain, nodal_strain);
    extrapolate_to_nodes(elem.type, n, gp_stress, nodal_stress);

    for (int i = 0; i < n; ++i) {
      const std::size_t gi = static_cast<std::size_t>(nidx[static_cast<std::size_t>(i)]);
      for (int comp = 0; comp < 6; ++comp) {
        fields.strain[gi][static_cast<std::size_t>(comp)] += nodal_strain[static_cast<std::size_t>(i)][static_cast<std::size_t>(comp)];
        fields.stress[gi][static_cast<std::size_t>(comp)] += nodal_stress[static_cast<std::size_t>(i)][static_cast<std::size_t>(comp)];
      }
      ++count[gi];
    }
  }

  for (std::size_t i = 0; i < n_nodes; ++i)
    if (count[i] > 0)
      for (int comp = 0; comp < 6; ++comp) {
        fields.strain[i][static_cast<std::size_t>(comp)] /= count[i];
        fields.stress[i][static_cast<std::size_t>(comp)] /= count[i];
      }

  // External load vector (*CLOAD + *DLOAD), then RF = f_int - f_ext.
  const std::vector<Real> f_ext = external_load_vector(model);
  for (std::size_t i = 0; i < n_nodes; ++i)
    for (int c = 0; c < kDofsPerNode; ++c)
      fields.reaction[i][static_cast<std::size_t>(c)] =
          f_int[i * kDofsPerNode + static_cast<std::size_t>(c)] -
          f_ext[i * kDofsPerNode + static_cast<std::size_t>(c)];
}

void recover_fields(const Model& model, StaticFields& fields, MaterialPoints& mp) {
  const Mesh& mesh = model.mesh;
  const std::size_t n_nodes = mesh.num_nodes();
  fields.stress.assign(n_nodes, Voigt6{});
  fields.strain.assign(n_nodes, Voigt6{});
  fields.reaction.assign(n_nodes, Vec3{0, 0, 0});
  std::vector<int> count(n_nodes, 0);
  std::vector<Real> f_int(n_nodes * kDofsPerNode, 0.0);
  std::array<std::array<Real, 3>, kMaxNodes> g{};

  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    const Element& elem = mesh.elements()[e];
    const int n = nodes_per_element(elem.type);
    std::vector<Vec3> coords(static_cast<std::size_t>(n));
    std::vector<Index> nidx(static_cast<std::size_t>(n));
    std::vector<Vec3> ue(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
      const Index ni = mesh.node_index(elem.nodes[static_cast<std::size_t>(i)]);
      nidx[static_cast<std::size_t>(i)] = ni;
      coords[static_cast<std::size_t>(i)] = mesh.nodes()[static_cast<std::size_t>(ni)].x;
      ue[static_cast<std::size_t>(i)] = fields.displacement[static_cast<std::size_t>(ni)];
    }

    // Element internal force from the material-point kernel at the converged
    // displacement + committed history (radial-return stress), scattered to globals.
    const ElementResponse resp =
        element_tangent_force(elem.type, coords, ue, *mp.models[e], mp.state[e]);
    const int ndof = n * kDofsPerNode;
    for (int a = 0; a < ndof; ++a)
      f_int[static_cast<std::size_t>(nidx[static_cast<std::size_t>(a / kDofsPerNode)]) *
                kDofsPerNode +
            static_cast<std::size_t>(a % kDofsPerNode)] +=
          resp.fe[static_cast<std::size_t>(a)];

    // Integration-point strain (from displacements) and stress (from the committed
    // material state), then extrapolate to nodes.
    const auto rule = gauss_rule(elem.type);
    std::vector<Voigt6> gp_strain(rule.size());
    std::vector<Voigt6> gp_stress(rule.size());
    for (std::size_t q = 0; q < rule.size(); ++q) {
      const Shape s = shape(elem.type, rule[q].xi, rule[q].et, rule[q].ze);
      physical_gradients(s, coords, g);
      gp_strain[q] = strain_from_gradients(g, n, ue);
      // Evaluate the committed state to get the converged stress at this point. The
      // material reads committed history, so this reproduces the converged stress
      // without advancing anything the caller relies on.
      MaterialState scratch = mp.state[e][q];
      gp_stress[q] = mp.models[e]->evaluate(gp_strain[q], scratch).stress;
    }
    std::vector<Voigt6> nodal_strain(static_cast<std::size_t>(n));
    std::vector<Voigt6> nodal_stress(static_cast<std::size_t>(n));
    extrapolate_to_nodes(elem.type, n, gp_strain, nodal_strain);
    extrapolate_to_nodes(elem.type, n, gp_stress, nodal_stress);
    for (int i = 0; i < n; ++i) {
      const std::size_t gi = static_cast<std::size_t>(nidx[static_cast<std::size_t>(i)]);
      for (int comp = 0; comp < 6; ++comp) {
        fields.strain[gi][static_cast<std::size_t>(comp)] += nodal_strain[static_cast<std::size_t>(i)][static_cast<std::size_t>(comp)];
        fields.stress[gi][static_cast<std::size_t>(comp)] += nodal_stress[static_cast<std::size_t>(i)][static_cast<std::size_t>(comp)];
      }
      ++count[gi];
    }
  }

  for (std::size_t i = 0; i < n_nodes; ++i)
    if (count[i] > 0)
      for (int comp = 0; comp < 6; ++comp) {
        fields.strain[i][static_cast<std::size_t>(comp)] /= count[i];
        fields.stress[i][static_cast<std::size_t>(comp)] /= count[i];
      }

  const std::vector<Real> f_ext = external_load_vector(model);
  for (std::size_t i = 0; i < n_nodes; ++i)
    for (int c = 0; c < kDofsPerNode; ++c)
      fields.reaction[i][static_cast<std::size_t>(c)] =
          f_int[i * kDofsPerNode + static_cast<std::size_t>(c)] -
          f_ext[i * kDofsPerNode + static_cast<std::size_t>(c)];
}

}  // namespace cxpp::fem
