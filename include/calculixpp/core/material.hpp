#pragma once
#include <optional>
#include <string>

#include "calculixpp/core/types.hpp"

namespace cxpp {

// Isotropic linear elasticity (spec: material-models — Phase 1).
struct ElasticIso {
  Real E{};   // Young's modulus
  Real nu{};  // Poisson's ratio
};

struct Material {
  std::string name;
  std::optional<ElasticIso> elastic;
  std::optional<Real> density;
};

}  // namespace cxpp
