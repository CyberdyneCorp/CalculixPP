#include "calculixpp/fem/stress.hpp"

#include <array>
#include <optional>
#include <span>
#include <vector>

#include "calculixpp/fem/element.hpp"
#include "calculixpp/fem/extrapolation.hpp"
#include "calculixpp/fem/loads.hpp"

namespace cxpp::fem {
namespace {

// Integration-point strain (engineering shear) and stress from displacements. The
// reported `strain` is the total mechanical strain B u; the stress uses the
// thermal-strain-corrected elastic strain sigma = D (eps_mech - eps_th), where
// `eps_th` is the (already-computed) thermal strain at this point (all-zero on the
// pure-mechanical path). (spec: heat-transfer — sigma = D (eps_mech - eps_th).)
void strain_stress_at(const std::array<std::array<Real, 3>, kMaxNodes>& g, int n,
                      const std::vector<Vec3>& ue, const D6& D, const Voigt6& eps_th,
                      Voigt6& strain, Voigt6& stress) {
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
    for (int j = 0; j < 6; ++j)
      v += D[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] *
           (strain[static_cast<std::size_t>(j)] - eps_th[static_cast<std::size_t>(j)]);
    stress[static_cast<std::size_t>(i)] = v;
  }
}

// Extrapolate integration-point Voigt6 values to element nodes via the shared scheme.
void extrapolate_to_nodes(ElementType type, int n, const std::vector<Voigt6>& gp,
                          std::vector<Voigt6>& nodal) {
  extrapolate_gp_to_nodes<6>(type, n, gp, nodal);
}

// Scatter one element's internal force f_int_e = Ke_e u_e - f_th_e into the global
// vector f_int. `f_th_e` is the element's thermal-strain load (empty span -> zero on
// the pure-mechanical path). Subtracting it makes f_int = ∫ Bᵀ D (eps_mech - eps_th),
// so the reaction RF = f_int - f_ext balances the thermal load correctly.
void accumulate_internal_force(ElementType type,
                               const std::vector<Vec3>& coords,
                               const std::vector<Index>& nidx,
                               const std::vector<Vec3>& ue, const ElasticIso& mat,
                               std::span<const Real> f_th_e,
                               std::vector<Real>& f_int) {
  const int n = nodes_per_element(type);
  const std::vector<Real> Ke = element_stiffness(type, coords, mat);
  const int ndof = n * kDofsPerNode;
  for (int a = 0; a < ndof; ++a) {
    Real v = 0.0;
    for (int b = 0; b < ndof; ++b)
      v += Ke[static_cast<std::size_t>(a) * static_cast<std::size_t>(ndof) + static_cast<std::size_t>(b)] *
           ue[static_cast<std::size_t>(b / kDofsPerNode)][static_cast<std::size_t>(b % kDofsPerNode)];
    if (!f_th_e.empty()) v -= f_th_e[static_cast<std::size_t>(a)];
    f_int[static_cast<std::size_t>(nidx[static_cast<std::size_t>(a / kDofsPerNode)]) * kDofsPerNode +
          static_cast<std::size_t>(a % kDofsPerNode)] += v;
  }
}

// Per-node temperature change (T - Tref) for element `e` given the model's applied
// temperature field and the element's optional *EXPANSION. Returns an empty vector
// when there is no expansion or no temperature change (so callers skip the thermal
// path). The values are NOT scaled by alpha (element_thermal_load applies alpha).
std::vector<Real> element_temp_change(const Model& model, const Element& elem,
                                      const std::optional<Expansion>& exp) {
  if (!exp || exp->empty() || model.applied_temperature.empty()) return {};
  const int n = nodes_per_element(elem.type);
  std::vector<Real> te(static_cast<std::size_t>(n), 0.0);
  bool any = false;
  for (int i = 0; i < n; ++i) {
    const auto it = model.applied_temperature.find(elem.nodes[static_cast<std::size_t>(i)]);
    if (it != model.applied_temperature.end()) {
      te[static_cast<std::size_t>(i)] = it->second - exp->t_ref;
      if (te[static_cast<std::size_t>(i)] != 0.0) any = true;
    }
  }
  return any ? te : std::vector<Real>{};
}

// Isotropic expansion coefficient of element `elem`, evaluated at its mean applied
// temperature for a temperature-dependent alpha(T) table (a constant table returns
// its single value regardless, so a constant deck is byte-for-byte unchanged). Nodes
// without an applied temperature contribute t_ref (their thermal strain is zero).
Real element_alpha(const Model& model, const Element& elem, const Expansion& exp) {
  const int n = nodes_per_element(elem.type);
  Real tsum = 0.0;
  for (int i = 0; i < n; ++i) {
    const auto it = model.applied_temperature.find(elem.nodes[static_cast<std::size_t>(i)]);
    tsum += (it == model.applied_temperature.end()) ? exp.t_ref : it->second;
  }
  return exp.alpha.at(n > 0 ? tsum / static_cast<Real>(n) : exp.t_ref);
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
  const std::vector<std::optional<Expansion>> expansion = model.element_expansion();
  const std::vector<bool> active = model.element_active_mask();
  std::vector<Vec3> coords, ue;
  std::vector<Index> nidx;
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    if (!active[e]) continue;  // *MODEL CHANGE, REMOVE: no internal force
    const Element& elem = mesh.elements()[e];
    gather_element(mesh, elem, u, coords, nidx, ue);
    const std::vector<Real> te = element_temp_change(model, elem, expansion[e]);
    std::vector<Real> f_th_e;
    if (!te.empty())
      f_th_e = element_thermal_load(elem.type, coords, elastic_iso_D(elastic[e]),
                                    element_alpha(model, elem, *expansion[e]), te);
    accumulate_internal_force(elem.type, coords, nidx, ue, elastic[e], f_th_e, f_int);
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
  const std::vector<std::optional<Expansion>> expansion = model.element_expansion();
  const std::vector<bool> active = model.element_active_mask();
  std::array<std::array<Real, 3>, kMaxNodes> g{};

  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    if (!active[e]) continue;  // *MODEL CHANGE, REMOVE: no stress / reaction contribution
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

    // Thermal strain: per-node temperature change (T - Tref) and its equivalent load.
    const std::vector<Real> te = element_temp_change(model, elem, expansion[e]);
    const Real alpha = te.empty() ? 0.0 : element_alpha(model, elem, *expansion[e]);
    std::vector<Real> f_th_e;
    if (!te.empty())
      f_th_e = element_thermal_load(elem.type, coords, D, alpha, te);

    // Internal force f_int_e = Ke_e u_e - f_th_e, scattered to global DOFs.
    accumulate_internal_force(elem.type, coords, nidx, ue, elastic[e], f_th_e, f_int);

    // Integration-point strains and stresses (stress uses eps_mech - eps_th).
    const auto rule = gauss_rule(elem.type);
    std::vector<Voigt6> gp_strain(rule.size());
    std::vector<Voigt6> gp_stress(rule.size());
    for (std::size_t q = 0; q < rule.size(); ++q) {
      const Shape s = shape(elem.type, rule[q].xi, rule[q].et, rule[q].ze);
      physical_gradients(s, coords, g);
      Real dT = 0.0;
      for (int k = 0; k < n; ++k)
        dT += s.N[static_cast<std::size_t>(k)] * (te.empty() ? 0.0 : te[static_cast<std::size_t>(k)]);
      const Real eth = alpha * dT;
      const Voigt6 eps_th{eth, eth, eth, 0.0, 0.0, 0.0};
      strain_stress_at(g, n, ue, D, eps_th, gp_strain[q], gp_stress[q]);
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
  const std::vector<bool> active = model.element_active_mask();
  std::array<std::array<Real, 3>, kMaxNodes> g{};

  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    if (!active[e]) continue;  // *MODEL CHANGE, REMOVE: no stress / reaction contribution
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
