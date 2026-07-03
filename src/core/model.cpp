#include "calculixpp/core/model.hpp"

#include <cctype>
#include <span>
#include <stdexcept>
#include <string>

#include "calculixpp/fem/element.hpp"

namespace cxpp {

std::optional<ElasticIso> Material::effective_elastic() const {
  if (elastic) return *elastic;
  if (hyperelastic && !hyperelastic->empty()) {
    // Derive (E, nu) from the neo-Hookean initial moduli mu = 2 C10, kappa = 2/D1:
    //   E = 9 K mu / (3 K + mu),  nu = (3 K - 2 mu) / (2 (3 K + mu)).
    const Real mu = hyperelastic->mu != 0.0 ? hyperelastic->mu
                                            : 2.0 * hyperelastic->c10;
    const Real K = hyperelastic->kappa != 0.0
                       ? hyperelastic->kappa
                       : (hyperelastic->d1 != 0.0 ? 2.0 / hyperelastic->d1 : 0.0);
    if (mu > 0.0 && K > 0.0) {
      const Real E = 9.0 * K * mu / (3.0 * K + mu);
      const Real nu = (3.0 * K - 2.0 * mu) / (2.0 * (3.0 * K + mu));
      return ElasticIso{E, nu};
    }
  }
  if (user && !user->empty() && user->constants.size() >= 2)
    return ElasticIso{user->constants[0], user->constants[1]};
  return std::nullopt;
}

std::vector<ElasticIso> Model::element_elastic() const {
  std::vector<ElasticIso> out(mesh.num_elements());
  std::vector<bool> assigned(mesh.num_elements(), false);

  for (const auto& section : sections) {
    const std::vector<Index>* ids = mesh.elset(section.elset);
    if (ids == nullptr) {
      throw std::runtime_error("solid section references unknown elset '" +
                               section.elset + "'");
    }
    const auto mat = materials.find(section.material);
    if (mat == materials.end()) {
      throw std::runtime_error("solid section references unknown material '" +
                               section.material + "'");
    }
    const std::optional<ElasticIso> eff = mat->second.effective_elastic();
    if (!eff) {
      throw std::runtime_error("material '" + section.material +
                               "' has no *ELASTIC (or hyperelastic/user) data");
    }
    for (const Index elem_id : *ids) {
      const Index ei = mesh.element_index(elem_id);
      if (ei < 0) {
        throw std::runtime_error("elset '" + section.elset +
                                 "' references unknown element");
      }
      out[static_cast<std::size_t>(ei)] = *eff;
      assigned[static_cast<std::size_t>(ei)] = true;
    }
  }

  for (std::size_t i = 0; i < assigned.size(); ++i) {
    if (!assigned[i]) {
      throw std::runtime_error("element " +
                               std::to_string(mesh.elements()[i].id) +
                               " has no solid section / material");
    }
  }
  return out;
}

std::vector<std::optional<Plastic>> Model::element_plastic() const {
  std::vector<std::optional<Plastic>> out(mesh.num_elements());
  for (const auto& section : sections) {
    const std::vector<Index>* ids = mesh.elset(section.elset);
    if (ids == nullptr) continue;  // element_elastic() already validates elsets
    const auto mat = materials.find(section.material);
    if (mat == materials.end()) continue;
    for (const Index elem_id : *ids) {
      const Index ei = mesh.element_index(elem_id);
      if (ei >= 0) out[static_cast<std::size_t>(ei)] = mat->second.plastic;
    }
  }
  return out;
}

namespace {
// Yield stress at accumulated equivalent plastic strain `ep`, linearly interpolated
// on the *PLASTIC hardening table (flat extrapolation past the ends). Mirrors the
// J2 material's `hardening()` so the deposited plastic work is consistent with the
// stress the mechanical solve return-mapped to.
Real yield_at(const Plastic& p, Real ep) {
  const std::size_t n = p.yield.size();
  if (n == 0) return 0.0;
  if (n == 1 || ep <= p.eqplastic[0]) return p.yield[0];
  if (ep >= p.eqplastic[n - 1]) return p.yield[n - 1];
  for (std::size_t i = 1; i < n; ++i)
    if (ep <= p.eqplastic[i]) {
      const Real t = (ep - p.eqplastic[i - 1]) /
                     (p.eqplastic[i] - p.eqplastic[i - 1]);
      return p.yield[i - 1] + t * (p.yield[i] - p.yield[i - 1]);
    }
  return p.yield[n - 1];
}
}  // namespace

std::vector<Real> Model::plastic_dissipation_heat(
    const std::vector<Real>& eqplastic_by_elem) const {
  std::vector<Real> heat(mesh.num_elements(), 0.0);
  if (taylor_quinney <= 0.0 || !has_plasticity()) return heat;
  const std::vector<std::optional<Plastic>> plas = element_plastic();
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    if (!plas[e] || plas[e]->empty()) continue;
    const Real ep = e < eqplastic_by_elem.size() ? eqplastic_by_elem[e] : 0.0;
    if (ep <= 0.0) continue;
    // Plastic work density w_p ~= sigma_y(ep) * ep (J2, associative flow). Deposit a
    // fraction beta of it over the element volume as heat.
    const Element& el = mesh.elements()[e];
    const int nn = nodes_per_element(el.type);
    std::vector<Vec3> coords(static_cast<std::size_t>(nn));
    for (int i = 0; i < nn; ++i) {
      const Index ni = mesh.node_index(el.nodes[static_cast<std::size_t>(i)]);
      coords[static_cast<std::size_t>(i)] = mesh.nodes()[static_cast<std::size_t>(ni)].x;
    }
    const Real vol = fem::element_volume(el.type, std::span<const Vec3>(coords));
    heat[e] = taylor_quinney * yield_at(*plas[e], ep) * ep * vol;
  }
  return heat;
}

std::vector<std::optional<Hyperelastic>> Model::element_hyperelastic() const {
  std::vector<std::optional<Hyperelastic>> out(mesh.num_elements());
  for (const auto& section : sections) {
    const std::vector<Index>* ids = mesh.elset(section.elset);
    if (ids == nullptr) continue;
    const auto mat = materials.find(section.material);
    if (mat == materials.end()) continue;
    for (const Index elem_id : *ids) {
      const Index ei = mesh.element_index(elem_id);
      if (ei >= 0) out[static_cast<std::size_t>(ei)] = mat->second.hyperelastic;
    }
  }
  return out;
}

std::vector<std::optional<UserMaterial>> Model::element_user_material() const {
  std::vector<std::optional<UserMaterial>> out(mesh.num_elements());
  for (const auto& section : sections) {
    const std::vector<Index>* ids = mesh.elset(section.elset);
    if (ids == nullptr) continue;
    const auto mat = materials.find(section.material);
    if (mat == materials.end()) continue;
    for (const Index elem_id : *ids) {
      const Index ei = mesh.element_index(elem_id);
      if (ei >= 0) out[static_cast<std::size_t>(ei)] = mat->second.user;
    }
  }
  return out;
}

std::vector<std::optional<Expansion>> Model::element_expansion() const {
  std::vector<std::optional<Expansion>> out(mesh.num_elements());
  for (const auto& section : sections) {
    const std::vector<Index>* ids = mesh.elset(section.elset);
    if (ids == nullptr) continue;  // element_elastic() already validates elsets
    const auto mat = materials.find(section.material);
    if (mat == materials.end()) continue;
    for (const Index elem_id : *ids) {
      const Index ei = mesh.element_index(elem_id);
      if (ei >= 0) out[static_cast<std::size_t>(ei)] = mat->second.expansion;
    }
  }
  return out;
}

bool Model::has_thermal_strain() const {
  if (applied_temperature.empty()) return false;
  for (const auto& [name, mat] : materials)
    if (mat.expansion && !mat.expansion->empty()) return true;
  return false;
}

bool Model::has_nonlinear_material() const {
  for (const auto& [name, mat] : materials) {
    if (mat.plastic && !mat.plastic->empty()) return true;
    if (mat.hyperelastic && !mat.hyperelastic->empty()) return true;
    if (mat.user && !mat.user->empty()) return true;
  }
  return false;
}

bool Model::has_plasticity() const {
  for (const auto& [name, mat] : materials)
    if (mat.plastic && !mat.plastic->empty()) return true;
  return false;
}

std::vector<bool> Model::element_active_mask() const {
  std::vector<bool> active(mesh.num_elements(), true);
  for (const Index elem_id : deactivated_elements) {
    const Index ei = mesh.element_index(elem_id);
    if (ei >= 0) active[static_cast<std::size_t>(ei)] = false;
  }
  return active;
}

SurfaceBehavior Model::contact_behavior(const ContactPair& pair) const {
  const auto it = surface_interactions.find(pair.interaction);
  if (it != surface_interactions.end() && it->second.has_behavior)
    return it->second.behavior;
  return SurfaceBehavior{};  // CalculiX default: hard pressure-overclosure
}

Friction Model::contact_friction(const ContactPair& pair) const {
  const auto it = surface_interactions.find(pair.interaction);
  if (it != surface_interactions.end() && it->second.friction.has)
    return it->second.friction;
  return Friction{};  // no *FRICTION -> frictionless normal contact
}

GapConductance Model::contact_conductance(const ContactPair& pair) const {
  const auto it = surface_interactions.find(pair.interaction);
  if (it != surface_interactions.end() && it->second.conductance.has)
    return it->second.conductance;
  return GapConductance{};  // no *GAP CONDUCTANCE -> thermally open interface
}

GapHeatGeneration Model::contact_heat_generation(const ContactPair& pair) const {
  const auto it = surface_interactions.find(pair.interaction);
  if (it != surface_interactions.end() && it->second.heat_generation.has)
    return it->second.heat_generation;
  return GapHeatGeneration{};  // no *GAP HEAT GENERATION -> no gap source
}

bool Model::has_thermal_contact() const {
  for (const ContactPair& pair : contact_pairs)
    if (contact_conductance(pair).has || contact_heat_generation(pair).has) return true;
  return false;
}

Model Model::with_active_contact_pairs() const {
  if (contact_pair_changes.empty()) return *this;  // no toggling -> every pair active
  // Case-insensitive surface-pair match, order-independent (the two data lines may list
  // the surfaces in either order).
  auto matches = [](const ContactPair& p, const ContactPairChange& c) {
    auto up = [](std::string s) {
      for (char& ch : s) ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
      return s;
    };
    const std::string ps = up(p.slave_surface), pm = up(p.master_surface);
    const std::string ca = up(c.surface_a), cb = up(c.surface_b);
    return (ps == ca && pm == cb) || (ps == cb && pm == ca);
  };
  Model m = *this;
  m.contact_pairs.clear();
  for (const ContactPair& pair : contact_pairs) {
    bool active = true;  // active unless the LAST matching change is a REMOVE
    for (const ContactPairChange& c : contact_pair_changes)
      if (matches(pair, c)) active = c.add;
    if (active) m.contact_pairs.push_back(pair);
  }
  return m;
}

std::vector<Real> Model::element_density() const {
  std::vector<Real> out(mesh.num_elements(), 0.0);
  for (const auto& section : sections) {
    const std::vector<Index>* ids = mesh.elset(section.elset);
    if (ids == nullptr) continue;  // element_elastic() already validates elsets
    const auto mat = materials.find(section.material);
    if (mat == materials.end() || !mat->second.density) continue;
    for (const Index elem_id : *ids) {
      const Index ei = mesh.element_index(elem_id);
      if (ei >= 0) out[static_cast<std::size_t>(ei)] = *mat->second.density;
    }
  }
  return out;
}

std::vector<Real> Model::element_conductivity() const {
  std::vector<Real> out(mesh.num_elements(), 0.0);
  std::vector<bool> assigned(mesh.num_elements(), false);

  for (const auto& section : sections) {
    const std::vector<Index>* ids = mesh.elset(section.elset);
    if (ids == nullptr) {
      throw std::runtime_error("solid section references unknown elset '" +
                               section.elset + "'");
    }
    const auto mat = materials.find(section.material);
    if (mat == materials.end()) {
      throw std::runtime_error("solid section references unknown material '" +
                               section.material + "'");
    }
    if (!mat->second.thermal || mat->second.thermal->conductivity.empty()) {
      throw std::runtime_error("material '" + section.material +
                               "' has no *CONDUCTIVITY data (needed for heat transfer)");
    }
    for (const Index elem_id : *ids) {
      const Index ei = mesh.element_index(elem_id);
      if (ei < 0) {
        throw std::runtime_error("elset '" + section.elset +
                                 "' references unknown element");
      }
      out[static_cast<std::size_t>(ei)] = mat->second.thermal->conductivity.first();
      assigned[static_cast<std::size_t>(ei)] = true;
    }
  }

  for (std::size_t i = 0; i < assigned.size(); ++i) {
    if (!assigned[i]) {
      throw std::runtime_error("element " +
                               std::to_string(mesh.elements()[i].id) +
                               " has no solid section / material");
    }
  }
  return out;
}

std::vector<Real> Model::element_heat_capacity() const {
  std::vector<Real> out(mesh.num_elements(), 0.0);
  for (const auto& section : sections) {
    const std::vector<Index>* ids = mesh.elset(section.elset);
    if (ids == nullptr) continue;  // element_conductivity() validates elsets
    const auto mat = materials.find(section.material);
    if (mat == materials.end() || !mat->second.thermal || !mat->second.density)
      continue;
    const Real rho_c = *mat->second.density * mat->second.thermal->specific_heat.first();
    for (const Index elem_id : *ids) {
      const Index ei = mesh.element_index(elem_id);
      if (ei >= 0) out[static_cast<std::size_t>(ei)] = rho_c;
    }
  }
  return out;
}

namespace {

// Mean of an element's nodal temperatures (used to evaluate temperature-dependent
// properties at a single per-element temperature — the conduction/capacitance
// kernels take one scalar property per element). `node_temp` is aligned with mesh
// node indices; missing entries fall back to 0.
Real element_mean_temp(const Mesh& mesh, const Element& elem,
                       const std::vector<Real>& node_temp) {
  const int n = nodes_per_element(elem.type);
  Real sum = 0.0;
  for (int i = 0; i < n; ++i) {
    const Index ni = mesh.node_index(elem.nodes[static_cast<std::size_t>(i)]);
    if (ni >= 0 && static_cast<std::size_t>(ni) < node_temp.size())
      sum += node_temp[static_cast<std::size_t>(ni)];
  }
  return n > 0 ? sum / static_cast<Real>(n) : 0.0;
}

// Resolve the (material) owning each element via the solid sections; invokes `fn`
// with (element index, material) for every element covered by a section whose
// material exists. Last section wins, matching element_elastic().
template <typename Fn>
void for_each_element_material(
    const Mesh& mesh,
    const std::unordered_map<std::string, Material>& materials,
    const std::vector<SolidSection>& sections, Fn&& fn) {
  for (const auto& section : sections) {
    const std::vector<Index>* ids = mesh.elset(section.elset);
    if (ids == nullptr) continue;
    const auto mat = materials.find(section.material);
    if (mat == materials.end()) continue;
    for (const Index elem_id : *ids) {
      const Index ei = mesh.element_index(elem_id);
      if (ei >= 0) fn(static_cast<std::size_t>(ei), mat->second);
    }
  }
}

}  // namespace

bool Model::has_temp_dependent_thermal() const {
  for (const auto& [name, mat] : materials) {
    if (mat.thermal &&
        (mat.thermal->conductivity.value.size() > 1 ||
         mat.thermal->specific_heat.value.size() > 1))
      return true;
    if (mat.expansion && mat.expansion->alpha.value.size() > 1) return true;
  }
  return false;
}

std::vector<Real> Model::element_conductivity_at(
    const std::vector<Real>& node_temp) const {
  std::vector<Real> out = element_conductivity();  // constant / first-value baseline
  for_each_element_material(
      mesh, materials, sections, [&](std::size_t ei, const Material& mat) {
        if (mat.thermal && mat.thermal->conductivity.value.size() > 1)
          out[ei] = mat.thermal->conductivity.at(
              element_mean_temp(mesh, mesh.elements()[ei], node_temp));
      });
  return out;
}

std::vector<Real> Model::element_heat_capacity_at(
    const std::vector<Real>& node_temp) const {
  std::vector<Real> out = element_heat_capacity();  // constant / first-value baseline
  for_each_element_material(
      mesh, materials, sections, [&](std::size_t ei, const Material& mat) {
        if (mat.thermal && mat.density &&
            mat.thermal->specific_heat.value.size() > 1)
          out[ei] = *mat.density *
                    mat.thermal->specific_heat.at(
                        element_mean_temp(mesh, mesh.elements()[ei], node_temp));
      });
  return out;
}

Real Model::amplitude_factor(const std::string& name, Real lambda) const {
  if (name.empty()) return lambda;  // default linear 0->1 ramp over the step
  const auto it = amplitudes.find(name);
  if (it == amplitudes.end()) return lambda;
  return it->second.value_at(lambda * increment.total);  // amplitude in step time
}

}  // namespace cxpp
