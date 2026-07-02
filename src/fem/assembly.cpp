#include "calculixpp/fem/assembly.hpp"

#include <cstdint>
#include <stdexcept>
#include <unordered_map>

#include "calculixpp/fem/element.hpp"
#include "calculixpp/fem/loads.hpp"

namespace cxpp::fem {

LinearSystem assemble_linear_static(const Model& model) {
  const Mesh& mesh = model.mesh;
  const auto n_nodes = mesh.num_nodes();
  const std::size_t n_dof = n_nodes * kDofsPerNode;

  LinearSystem sys;
  sys.dof_eq.assign(n_dof, -1);
  sys.prescribed.assign(n_dof, 0.0);
  std::vector<bool> constrained(n_dof, false);

  // Apply single-point constraints.
  for (const Spc& spc : model.spcs) {
    const Index ni = mesh.node_index(spc.node_id);
    if (ni < 0) throw std::runtime_error("*BOUNDARY references unknown node");
    if (spc.comp < 1 || spc.comp > kDofsPerNode)
      throw std::runtime_error("*BOUNDARY dof out of range");
    const std::size_t g = static_cast<std::size_t>(ni) * kDofsPerNode +
                          static_cast<std::size_t>(spc.comp - 1);
    constrained[g] = true;
    sys.prescribed[g] = spc.value;
  }

  // Number the free DOFs.
  Index n_free = 0;
  for (std::size_t g = 0; g < n_dof; ++g) {
    if (!constrained[g]) sys.dof_eq[g] = n_free++;
  }
  sys.n_free = n_free;
  sys.rhs.assign(static_cast<std::size_t>(n_free), 0.0);

  const std::vector<ElasticIso> elastic = model.element_elastic();

  // Accumulate the free/free stiffness (duplicates summed) and the eliminated RHS.
  std::unordered_map<std::int64_t, Real> kmap;
  const auto key = [n_free](Index r, Index c) {
    return static_cast<std::int64_t>(r) * n_free + c;
  };

  std::vector<Vec3> coords;
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    const Element& elem = mesh.elements()[e];
    const int n = nodes_per_element(elem.type);
    coords.resize(static_cast<std::size_t>(n));
    std::vector<Index> node_idx(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
      const Index ni = mesh.node_index(elem.nodes[static_cast<std::size_t>(i)]);
      if (ni < 0) throw std::runtime_error("element references unknown node");
      node_idx[static_cast<std::size_t>(i)] = ni;
      coords[static_cast<std::size_t>(i)] = mesh.nodes()[static_cast<std::size_t>(ni)].x;
    }

    const std::vector<Real> Ke =
        element_stiffness(elem.type, coords, elastic[e]);
    const int ndof = n * kDofsPerNode;

    for (int a = 0; a < ndof; ++a) {
      const std::size_t gA = static_cast<std::size_t>(node_idx[static_cast<std::size_t>(a / kDofsPerNode)]) *
                                 kDofsPerNode + static_cast<std::size_t>(a % kDofsPerNode);
      const Index eqA = sys.dof_eq[gA];
      if (eqA < 0) continue;  // constrained row: no equation
      for (int b = 0; b < ndof; ++b) {
        const std::size_t gB = static_cast<std::size_t>(node_idx[static_cast<std::size_t>(b / kDofsPerNode)]) *
                                   kDofsPerNode + static_cast<std::size_t>(b % kDofsPerNode);
        const Real kab = Ke[static_cast<std::size_t>(a) * static_cast<std::size_t>(ndof) + static_cast<std::size_t>(b)];
        const Index eqB = sys.dof_eq[gB];
        if (eqB >= 0) {
          kmap[key(eqA, eqB)] += kab;
        } else {
          // constrained column -> move to RHS: -K_fc u_c
          sys.rhs[static_cast<std::size_t>(eqA)] -= kab * sys.prescribed[gB];
        }
      }
    }
  }

  // External loads (*CLOAD + *DLOAD pressure) on free DOFs.
  const std::vector<Real> fext = external_load_vector(model);
  for (std::size_t g = 0; g < n_dof; ++g) {
    const Index eq = sys.dof_eq[g];
    if (eq >= 0) sys.rhs[static_cast<std::size_t>(eq)] += fext[g];
  }

  sys.rows.reserve(kmap.size());
  sys.cols.reserve(kmap.size());
  sys.vals.reserve(kmap.size());
  for (const auto& [k, v] : kmap) {
    sys.rows.push_back(static_cast<Index>(k / n_free));
    sys.cols.push_back(static_cast<Index>(k % n_free));
    sys.vals.push_back(v);
  }
  return sys;
}

}  // namespace cxpp::fem
