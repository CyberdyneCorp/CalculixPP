#include "calculixpp/fem/loads.hpp"

#include <array>
#include <stdexcept>

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

}  // namespace

std::vector<int> face_nodes(ElementType type, int face) {
  if (face < 1 || face > 4) throw std::runtime_error("tet face must be 1..4");
  const int nf = (type == ElementType::C3D4) ? 3 : 6;
  std::vector<int> out(static_cast<std::size_t>(nf));
  for (int i = 0; i < nf; ++i) out[static_cast<std::size_t>(i)] = kFace[face - 1][i];
  return out;
}

std::vector<Real> external_load_vector(const Model& model) {
  const Mesh& mesh = model.mesh;
  std::vector<Real> f(mesh.num_nodes() * kDofsPerNode, 0.0);

  for (const Cload& cl : model.cloads) {
    const Index ni = mesh.node_index(cl.node_id);
    if (ni >= 0 && cl.comp >= 1 && cl.comp <= kDofsPerNode)
      f[static_cast<std::size_t>(ni) * kDofsPerNode + static_cast<std::size_t>(cl.comp - 1)] += cl.value;
  }

  for (const Dload& dl : model.dloads) {
    const Index ei = mesh.element_index(dl.elem_id);
    if (ei < 0) throw std::runtime_error("*DLOAD references unknown element");
    const Element& elem = mesh.elements()[static_cast<std::size_t>(ei)];
    const std::vector<int> fn = face_nodes(elem.type, dl.face);
    const int nf = static_cast<int>(fn.size());

    // Physical coordinates of the face nodes.
    std::vector<Vec3> x(static_cast<std::size_t>(nf));
    std::vector<Index> gnode(static_cast<std::size_t>(nf));
    for (int i = 0; i < nf; ++i) {
      const Index ni = mesh.node_index(elem.nodes[static_cast<std::size_t>(fn[static_cast<std::size_t>(i)])]);
      gnode[static_cast<std::size_t>(i)] = ni;
      x[static_cast<std::size_t>(i)] = mesh.nodes()[static_cast<std::size_t>(ni)].x;
    }

    // F_i = -p * integral( N_i * (t_xi x t_eta) ) over the face.
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
              -dl.pressure * N[static_cast<std::size_t>(i)] * nrm[static_cast<std::size_t>(d)] * gp.w;
    }

    for (int i = 0; i < nf; ++i)
      for (int d = 0; d < 3; ++d)
        f[static_cast<std::size_t>(gnode[static_cast<std::size_t>(i)]) * kDofsPerNode + static_cast<std::size_t>(d)] +=
            force[static_cast<std::size_t>(i)][static_cast<std::size_t>(d)];
  }
  return f;
}

}  // namespace cxpp::fem
