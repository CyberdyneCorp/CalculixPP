#pragma once
#include <cstddef>
#include <vector>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/types.hpp"

// Stress-based high-cycle-fatigue (HCF) evaluation (spec: high-cycle-fatigue — *HCF).
// Treats the recovered stress field of a preceding stress-producing model as the cyclic
// stress amplitude at each node, reduces it to a scalar uniaxial-equivalent amplitude S_a
// (the model's FatigueCriterion), and inverts the material Basquin S-N curve
//   S_a = a * N^b   ->   N = (S_a / a)^(1/b)
// for the cycles-to-failure at every node. The worst-case (critical) location is the node
// of largest S_a, hence smallest life. This is a stress-life estimate — NOT the CalculiX
// crack-growth cumulative HCF path (hcfs.f/combilcfhcf.f).
namespace cxpp::numerics {

// Result of a *HCF evaluation. Per-node arrays are aligned with the mesh node indices;
// `life[i]` is the cycles-to-failure and `amplitude[i]` the scalar stress amplitude at
// node i. The worst-case block names the critical node (smallest life).
struct HcfReport {
  std::vector<Real> life;       // cycles-to-failure N per node
  std::vector<Real> amplitude;  // scalar stress amplitude S_a per node

  Index worst_node_id{};    // user id of the critical (smallest-life) node
  std::size_t worst_index{0};  // its 0-based mesh node index
  Vec3 worst_location{0, 0, 0};  // its coordinates
  Real worst_amplitude{0.0};   // its stress amplitude S_a
  Real worst_life{0.0};        // its cycles-to-failure N
};

// Run the stress-life HCF evaluation for `model`. Recovers the stress field with the
// linear-static path, applies the model's FatigueCriterion + the (single) material S-N
// curve, and returns the per-node life/amplitude + worst-case location.
//   Throws std::runtime_error when no material carries an S-N (*FATIGUE) curve, or when the
//   model has no elastic material to recover a stress field from (missing source).
HcfReport evaluate_hcf(const Model& model);

}  // namespace cxpp::numerics
