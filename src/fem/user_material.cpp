// User-material plug-in registry and a reference user law (spec: material-models 4.6 —
// *USER MATERIAL C++ interface, *DEPVAR state, *RATEDEPENDENT scaling).
//
// The registry is the extensibility seam: a *USER MATERIAL card names a factory
// registered here, and make_material_points (fem/assembly.cpp) builds the model from
// the parsed UserMaterial. LinearUserMaterial is a self-contained reference law that
// proves the whole path end to end without shipping a third-party constitutive model:
// it is a linear-elastic law parameterized purely by the *USER MATERIAL constants
// {E, nu}, records a *DEPVAR state variable, and honors a *RATEDEPENDENT scale.
#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>

#include "calculixpp/fem/element.hpp"
#include "calculixpp/fem/material_model.hpp"

namespace cxpp::fem {
namespace {

std::string upper_key(const std::string& s) {
  std::string out = s;
  std::transform(out.begin(), out.end(), out.begin(),
                 [](unsigned char c) { return std::toupper(c); });
  return out;
}

}  // namespace

UserMaterialRegistry& UserMaterialRegistry::instance() {
  static UserMaterialRegistry reg;
  return reg;
}

bool UserMaterialRegistry::register_factory(const std::string& name,
                                            UserMaterialFactory factory) {
  factories_[upper_key(name)] = std::move(factory);
  return true;
}

const UserMaterialFactory* UserMaterialRegistry::find(
    const std::string& name) const {
  const auto it = factories_.find(upper_key(name));
  return it == factories_.end() ? nullptr : &it->second;
}

bool register_user_material(const std::string& name, UserMaterialFactory factory) {
  return UserMaterialRegistry::instance().register_factory(name, std::move(factory));
}

// ---------------------------------------------------------------------------------
// Reference user law: linear isotropic elasticity driven by *USER MATERIAL constants.
LinearUserMaterial::LinearUserMaterial(Real E, Real nu, Real rate_scale)
    : D_(elastic_iso_D(ElasticIso{E, nu})), rate_scale_(rate_scale) {
  if (rate_scale_ == 0.0) rate_scale_ = 1.0;
}

MaterialResponse LinearUserMaterial::evaluate(const Voigt6& strain,
                                              MaterialState& state) const {
  MaterialResponse r;
  for (int i = 0; i < kVoigt6; ++i) {
    Real v = 0.0;
    for (int j = 0; j < kVoigt6; ++j)
      v += D_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] *
           strain[static_cast<std::size_t>(j)];
    r.stress[static_cast<std::size_t>(i)] = rate_scale_ * v;
    for (int j = 0; j < kVoigt6; ++j)
      r.tangent[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
          rate_scale_ * D_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
  }
  // Exercise the *DEPVAR slot: store a scalar state variable (here the equivalent
  // strain) so the trial->committed promotion path is proven for user models too.
  if (!state.depvar.empty()) {
    const Real eq = std::sqrt(strain[0] * strain[0] + strain[1] * strain[1] +
                              strain[2] * strain[2]);
    state.depvar[0] = eq;
  }
  return r;
}

bool LinearUserMaterial::register_default(const std::string& name) {
  return register_user_material(
      name, [](const UserMaterial& um, const ElasticIso& el) {
        // Constants convention: {E, nu}. Fall back to the material's *ELASTIC block
        // when *USER MATERIAL provided no constants (lets a deck combine both).
        Real E = el.E, nu = el.nu;
        if (um.constants.size() >= 1) E = um.constants[0];
        if (um.constants.size() >= 2) nu = um.constants[1];
        if (E <= 0.0)
          throw std::runtime_error(
              "LINEAR_ELASTIC_USER needs E>0 from *USER MATERIAL constants or *ELASTIC");
        const Real scale = um.rate_scale.value_or(1.0);
        return std::make_unique<LinearUserMaterial>(E, nu, scale);
      });
}

namespace {
// Auto-register the reference law at load time so a deck can use it out of the box.
const bool kRegisteredDefaultUserMaterial = LinearUserMaterial::register_default();
}  // namespace

}  // namespace cxpp::fem
