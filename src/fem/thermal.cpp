#include "calculixpp/fem/thermal.hpp"

#include <cstdint>
#include <stdexcept>
#include <unordered_map>

#include "calculixpp/fem/element.hpp"
#include "calculixpp/fem/loads.hpp"

namespace cxpp::fem {
namespace {

// Scatter a dense n x n element block into the reduced free/free COO map, moving
// prescribed-temperature columns to the rhs (q_f -= Kt_fc * T_c). Scalar analogue
// of the mechanical scatter_block, without a constraint transform. Shared by the
// conduction and film blocks.
void scatter_block(const std::vector<Real>& block, int n,
                   const std::vector<Index>& nidx, const LinearSystem& sys,
                   Index n_free, std::unordered_map<std::int64_t, Real>& kmap,
                   std::vector<Real>& rhs) {
  const auto key = [n_free](Index r, Index c) {
    return static_cast<std::int64_t>(r) * n_free + c;
  };
  for (int a = 0; a < n; ++a) {
    const Index eqa = sys.dof_eq[static_cast<std::size_t>(nidx[static_cast<std::size_t>(a)])];
    if (eqa < 0) continue;  // prescribed row -> its equation is eliminated
    for (int b = 0; b < n; ++b) {
      const Real kab = block[static_cast<std::size_t>(a) * static_cast<std::size_t>(n) +
                             static_cast<std::size_t>(b)];
      if (kab == 0.0) continue;
      const std::size_t nb = static_cast<std::size_t>(nidx[static_cast<std::size_t>(b)]);
      const Index eqb = sys.dof_eq[nb];
      if (eqb >= 0) {
        kmap[key(eqa, eqb)] += kab;
      } else {
        rhs[static_cast<std::size_t>(eqa)] -= kab * sys.prescribed[nb];
      }
    }
  }
}

// Flush the accumulated free/free COO map into the LinearSystem triplet arrays.
void flush_coo(const std::unordered_map<std::int64_t, Real>& kmap, Index n_free,
               LinearSystem& sys) {
  sys.rows.reserve(kmap.size());
  sys.cols.reserve(kmap.size());
  sys.vals.reserve(kmap.size());
  for (const auto& [k, v] : kmap) {
    sys.rows.push_back(static_cast<Index>(k / n_free));
    sys.cols.push_back(static_cast<Index>(k % n_free));
    sys.vals.push_back(v);
  }
}

// Gather one element's node indices and coordinates.
void gather_element(const Mesh& mesh, const Element& elem, int n,
                    std::vector<Index>& nidx, std::vector<Vec3>& coords) {
  nidx.resize(static_cast<std::size_t>(n));
  coords.resize(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    const Index ni = mesh.node_index(elem.nodes[static_cast<std::size_t>(i)]);
    if (ni < 0) throw std::runtime_error("element references unknown node");
    nidx[static_cast<std::size_t>(i)] = ni;
    coords[static_cast<std::size_t>(i)] = mesh.nodes()[static_cast<std::size_t>(ni)].x;
  }
}

}  // namespace

void build_thermal_dof_map(const Model& model, LinearSystem& sys) {
  const Mesh& mesh = model.mesh;
  const std::size_t n_nodes = mesh.num_nodes();
  sys.prescribed.assign(n_nodes, 0.0);
  std::vector<bool> fixed(n_nodes, false);

  for (const TempBc& bc : model.temp_bcs) {
    const Index ni = mesh.node_index(bc.node_id);
    if (ni < 0) throw std::runtime_error("thermal *BOUNDARY references unknown node");
    fixed[static_cast<std::size_t>(ni)] = true;
    sys.prescribed[static_cast<std::size_t>(ni)] = bc.value;  // last writer wins
  }

  sys.dof_eq.assign(n_nodes, -1);
  Index eq = 0;
  for (std::size_t ni = 0; ni < n_nodes; ++ni)
    if (!fixed[ni]) sys.dof_eq[ni] = eq++;
  sys.n_free = eq;
  sys.rhs.assign(static_cast<std::size_t>(sys.n_free), 0.0);
}

LinearSystem assemble_conduction(const Model& model) {
  const Mesh& mesh = model.mesh;
  LinearSystem sys;
  build_thermal_dof_map(model, sys);
  const Index n_free = sys.n_free;

  const std::vector<Real> cond = model.element_conductivity();
  const std::vector<bool> active = model.element_active_mask();
  std::unordered_map<std::int64_t, Real> kmap;
  std::vector<Index> nidx;
  std::vector<Vec3> coords;

  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    if (!active[e]) continue;  // *MODEL CHANGE, REMOVE: element carries no conduction
    const Element& elem = mesh.elements()[e];
    const int n = nodes_per_element(elem.type);
    gather_element(mesh, elem, n, nidx, coords);
    const std::vector<Real> Kt = element_conduction(elem.type, coords, cond[e]);
    scatter_block(Kt, n, nidx, sys, n_free, kmap, sys.rhs);
  }

  // Film (*FILM): each face adds h ∫N_a N_b to the matrix and h T_sink ∫N_a to rhs.
  std::vector<Index> fnodes;
  for (const Film& fm : model.films) {
    const Index ei = mesh.element_index(fm.elem_id);
    if (ei < 0) throw std::runtime_error("*FILM references unknown element");
    const FaceSurface fs = face_surface_integrals(mesh, ei, fm.face);
    const int nf = static_cast<int>(fs.gnode.size());
    std::vector<Real> hM(fs.mass.size());
    for (std::size_t i = 0; i < hM.size(); ++i) hM[i] = fm.h * fs.mass[i];
    scatter_block(hM, nf, fs.gnode, sys, n_free, kmap, sys.rhs);
    for (int a = 0; a < nf; ++a) {
      const Index eqa = sys.dof_eq[static_cast<std::size_t>(fs.gnode[static_cast<std::size_t>(a)])];
      if (eqa >= 0)
        sys.rhs[static_cast<std::size_t>(eqa)] +=
            fm.h * fm.sink_temp * fs.load[static_cast<std::size_t>(a)];
    }
  }

  // Concentrated + distributed flux into the free rhs. Fluxes on prescribed nodes
  // become reactions (recovered as RFL at solve time), so they are dropped here.
  const std::vector<Real> q = thermal_load_vector(model);
  for (std::size_t ni = 0; ni < mesh.num_nodes(); ++ni) {
    const Index eq = sys.dof_eq[ni];
    if (eq >= 0) sys.rhs[static_cast<std::size_t>(eq)] += q[ni];
  }

  flush_coo(kmap, n_free, sys);
  return sys;
}

FullThermalSystem assemble_full_thermal(const Model& model, Real lambda) {
  const Mesh& mesh = model.mesh;
  const std::vector<Real> cond = model.element_conductivity();
  const std::vector<Real> rho_c = model.element_heat_capacity();
  const std::vector<bool> active = model.element_active_mask();
  FullThermalSystem out;

  std::vector<Index> nidx;
  std::vector<Vec3> coords;
  const auto emit = [](std::vector<Index>& rows, std::vector<Index>& cols,
                       std::vector<Real>& vals, Index r, Index c, Real v) {
    if (v == 0.0) return;
    rows.push_back(r);
    cols.push_back(c);
    vals.push_back(v);
  };

  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    if (!active[e]) continue;  // *MODEL CHANGE, REMOVE: no conduction / capacitance
    const Element& elem = mesh.elements()[e];
    const int n = nodes_per_element(elem.type);
    gather_element(mesh, elem, n, nidx, coords);
    const std::vector<Real> Kt = element_conduction(elem.type, coords, cond[e]);
    const std::vector<Real> Ce = element_capacitance(elem.type, coords, rho_c[e]);
    for (int a = 0; a < n; ++a)
      for (int b = 0; b < n; ++b) {
        const std::size_t idx = static_cast<std::size_t>(a) * static_cast<std::size_t>(n) +
                                static_cast<std::size_t>(b);
        emit(out.k_rows, out.k_cols, out.k_vals, nidx[static_cast<std::size_t>(a)],
             nidx[static_cast<std::size_t>(b)], Kt[idx]);
        emit(out.c_rows, out.c_cols, out.c_vals, nidx[static_cast<std::size_t>(a)],
             nidx[static_cast<std::size_t>(b)], Ce[idx]);
      }
  }

  out.flux = thermal_load_vector(model, lambda);

  // Film folds into the conductivity operator and its source into flux.
  for (const Film& fm : model.films) {
    const Index ei = mesh.element_index(fm.elem_id);
    if (ei < 0) throw std::runtime_error("*FILM references unknown element");
    const FaceSurface fs = face_surface_integrals(mesh, ei, fm.face);
    const int nf = static_cast<int>(fs.gnode.size());
    const Real s = model.amplitude_factor(fm.amplitude, lambda);
    for (int a = 0; a < nf; ++a) {
      out.flux[static_cast<std::size_t>(fs.gnode[static_cast<std::size_t>(a)])] +=
          fm.h * (s * fm.sink_temp) * fs.load[static_cast<std::size_t>(a)];
      for (int b = 0; b < nf; ++b)
        emit(out.k_rows, out.k_cols, out.k_vals, fs.gnode[static_cast<std::size_t>(a)],
             fs.gnode[static_cast<std::size_t>(b)],
             fm.h * fs.mass[static_cast<std::size_t>(a) * static_cast<std::size_t>(nf) +
                            static_cast<std::size_t>(b)]);
    }
  }
  return out;
}

}  // namespace cxpp::fem
