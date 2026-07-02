#pragma once
#include <array>
#include <vector>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/types.hpp"

namespace cxpp::fem {

// Local element node indices (0-based) of a tetrahedral face (1..4), CalculiX
// ifacet ordering: 3 corners for C3D4, 3 corners + 3 midsides for C3D10.
std::vector<int> face_nodes(ElementType type, int face);

// Global external load vector (size 3*num_nodes): concentrated loads (*CLOAD)
// plus consistent nodal forces from pressure (*DLOAD, P<face>).
std::vector<Real> external_load_vector(const Model& model);

}  // namespace cxpp::fem
