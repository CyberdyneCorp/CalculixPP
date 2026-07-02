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

// Distributed pressure on an element face (*DLOAD, P<face>). Positive pressure
// acts into the element (opposite the outward face normal). Face is 1..4 (tet).
struct Dload {
  Index elem_id{};
  int face{};      // 1..4
  Real pressure{0.0};
};

// The assembled analysis model for the linear-static slice.
class Model {
 public:
  Mesh mesh;
  std::unordered_map<std::string, Material> materials;
  std::vector<SolidSection> sections;
  std::vector<Spc> spcs;
  std::vector<Cload> cloads;
  std::vector<Dload> dloads;

  // Elastic properties per element (aligned with mesh.elements()), resolved from
  // the solid sections. Throws std::runtime_error on a missing elset/material or an
  // element left without a section.
  std::vector<ElasticIso> element_elastic() const;
};

}  // namespace cxpp
