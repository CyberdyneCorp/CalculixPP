#include "calculixpp/core/model.hpp"

#include <stdexcept>
#include <string>

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
    if (!mat->second.thermal || mat->second.thermal->conductivity == 0.0) {
      throw std::runtime_error("material '" + section.material +
                               "' has no *CONDUCTIVITY data (needed for heat transfer)");
    }
    for (const Index elem_id : *ids) {
      const Index ei = mesh.element_index(elem_id);
      if (ei < 0) {
        throw std::runtime_error("elset '" + section.elset +
                                 "' references unknown element");
      }
      out[static_cast<std::size_t>(ei)] = mat->second.thermal->conductivity;
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
    const Real rho_c = *mat->second.density * mat->second.thermal->specific_heat;
    for (const Index elem_id : *ids) {
      const Index ei = mesh.element_index(elem_id);
      if (ei >= 0) out[static_cast<std::size_t>(ei)] = rho_c;
    }
  }
  return out;
}

Real Model::amplitude_factor(const std::string& name, Real lambda) const {
  if (name.empty()) return lambda;  // default linear 0->1 ramp over the step
  const auto it = amplitudes.find(name);
  if (it == amplitudes.end()) return lambda;
  return it->second.value_at(lambda * increment.total);  // amplitude in step time
}

}  // namespace cxpp
