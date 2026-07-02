#include "calculixpp/fem/stress.hpp"

#include <array>
#include <vector>

#include "calculixpp/fem/element.hpp"

namespace cxpp::fem {
namespace {

// C3D10 corner extrapolation matrix a4 (inverse of the linear-tet shape functions
// at the 4 Gauss points) and midside adjacency (nonei10), from extrapolate.f.
constexpr Real kA4diag = 1.927050983124842;
constexpr Real kA4off = -0.309016994374947;
constexpr int kMid[6][2] = {{0, 1}, {1, 2}, {2, 0}, {0, 3}, {1, 3}, {2, 3}};

Voigt6 stress_at(const std::array<std::array<Real, 3>, kMaxNodes>& g, int n,
                 const std::vector<Vec3>& ue, const D6& D) {
  // Displacement gradient dudx[i][j] = sum_k g_k[j] * u_k[i].
  std::array<std::array<Real, 3>, 3> du{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      Real v = 0.0;
      for (int k = 0; k < n; ++k) v += g[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)] * ue[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)];
      du[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] = v;
    }
  const std::array<Real, 6> e{du[0][0], du[1][1], du[2][2],
                              du[0][1] + du[1][0], du[0][2] + du[2][0],
                              du[1][2] + du[2][1]};
  Voigt6 s{};
  for (int i = 0; i < 6; ++i) {
    Real v = 0.0;
    for (int j = 0; j < 6; ++j) v += D[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] * e[static_cast<std::size_t>(j)];
    s[static_cast<std::size_t>(i)] = v;
  }
  return s;
}

}  // namespace

void recover_fields(const Model& model, StaticFields& fields) {
  const Mesh& mesh = model.mesh;
  const std::size_t n_nodes = mesh.num_nodes();
  fields.stress.assign(n_nodes, Voigt6{});
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
    const std::vector<Real> Ke = element_stiffness(elem.type, coords, elastic[e]);
    const int ndof = n * kDofsPerNode;
    for (int a = 0; a < ndof; ++a) {
      Real v = 0.0;
      for (int b = 0; b < ndof; ++b)
        v += Ke[static_cast<std::size_t>(a) * static_cast<std::size_t>(ndof) + static_cast<std::size_t>(b)] *
             ue[static_cast<std::size_t>(b / kDofsPerNode)][static_cast<std::size_t>(b % kDofsPerNode)];
      f_int[static_cast<std::size_t>(nidx[static_cast<std::size_t>(a / kDofsPerNode)]) * kDofsPerNode +
            static_cast<std::size_t>(a % kDofsPerNode)] += v;
    }

    // Integration-point stresses.
    const auto rule = gauss_rule(elem.type);
    std::vector<Voigt6> gp_stress(rule.size());
    for (std::size_t q = 0; q < rule.size(); ++q) {
      const Shape s = shape(elem.type, rule[q].xi, rule[q].et, rule[q].ze);
      physical_gradients(s, coords, g);
      gp_stress[q] = stress_at(g, n, ue, D);
    }

    // Extrapolate to element nodes.
    std::vector<Voigt6> nodal(static_cast<std::size_t>(n));
    if (elem.type == ElementType::C3D4) {
      for (int i = 0; i < 4; ++i) nodal[static_cast<std::size_t>(i)] = gp_stress[0];
    } else {
      for (int c = 0; c < 4; ++c)
        for (int comp = 0; comp < 6; ++comp) {
          Real v = 0.0;
          for (int l = 0; l < 4; ++l)
            v += (c == l ? kA4diag : kA4off) * gp_stress[static_cast<std::size_t>(l)][static_cast<std::size_t>(comp)];
          nodal[static_cast<std::size_t>(c)][static_cast<std::size_t>(comp)] = v;
        }
      for (int mnode = 0; mnode < 6; ++mnode)
        for (int comp = 0; comp < 6; ++comp)
          nodal[static_cast<std::size_t>(4 + mnode)][static_cast<std::size_t>(comp)] =
              0.5 * (nodal[static_cast<std::size_t>(kMid[mnode][0])][static_cast<std::size_t>(comp)] +
                     nodal[static_cast<std::size_t>(kMid[mnode][1])][static_cast<std::size_t>(comp)]);
    }

    for (int i = 0; i < n; ++i) {
      const std::size_t gi = static_cast<std::size_t>(nidx[static_cast<std::size_t>(i)]);
      for (int comp = 0; comp < 6; ++comp) fields.stress[gi][static_cast<std::size_t>(comp)] += nodal[static_cast<std::size_t>(i)][static_cast<std::size_t>(comp)];
      ++count[gi];
    }
  }

  for (std::size_t i = 0; i < n_nodes; ++i)
    if (count[i] > 0)
      for (int comp = 0; comp < 6; ++comp) fields.stress[i][static_cast<std::size_t>(comp)] /= count[i];

  // External load vector, then RF = f_int - f_ext.
  std::vector<Real> f_ext(n_nodes * kDofsPerNode, 0.0);
  for (const Cload& cl : model.cloads) {
    const Index ni = mesh.node_index(cl.node_id);
    if (ni >= 0 && cl.comp >= 1 && cl.comp <= kDofsPerNode)
      f_ext[static_cast<std::size_t>(ni) * kDofsPerNode + static_cast<std::size_t>(cl.comp - 1)] += cl.value;
  }
  for (std::size_t i = 0; i < n_nodes; ++i)
    for (int c = 0; c < kDofsPerNode; ++c)
      fields.reaction[i][static_cast<std::size_t>(c)] =
          f_int[i * kDofsPerNode + static_cast<std::size_t>(c)] -
          f_ext[i * kDofsPerNode + static_cast<std::size_t>(c)];
}

}  // namespace cxpp::fem
