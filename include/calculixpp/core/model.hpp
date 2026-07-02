#pragma once
#include <string>
#include <unordered_map>
#include <vector>

#include "calculixpp/core/material.hpp"
#include "calculixpp/core/mesh.hpp"
#include "calculixpp/core/section.hpp"

namespace cxpp {

// Single-point constraint: prescribed nodal DOF (spec: loads-and-boundary-conditions).
struct Spc {
  Index node_id{};
  int comp{};      // 1..3 (x,y,z)
  Real value{0.0};
};

// Concentrated load on a nodal DOF (*CLOAD).
struct Cload {
  Index node_id{};
  int comp{};      // 1..3
  Real value{0.0};
};

// The assembled analysis model for the linear-static slice.
class Model {
 public:
  Mesh mesh;
  std::unordered_map<std::string, Material> materials;
  std::vector<SolidSection> sections;
  std::vector<Spc> spcs;
  std::vector<Cload> cloads;

  // Elastic properties per element (aligned with mesh.elements()), resolved from
  // the solid sections. Throws std::runtime_error on a missing elset/material or an
  // element left without a section.
  std::vector<ElasticIso> element_elastic() const;
};

}  // namespace cxpp
