#pragma once
#include <string_view>
#include <vector>

#include "calculixpp/core/types.hpp"

namespace cxpp {

// Solid element library.
//   Phase 1: linear/quadratic tetrahedra (C3D4, C3D10).
//   Phase 2 (workstream 3): hexahedra and wedges with full integration
//   (C3D8, C3D20, C3D6, C3D15) and reduced-integration variants (C3D8R, C3D20R).
enum class ElementType {
  C3D4,
  C3D10,
  C3D8,
  C3D8R,
  C3D20,
  C3D20R,
  C3D6,
  C3D15,
};

constexpr int nodes_per_element(ElementType t) {
  switch (t) {
    case ElementType::C3D4:
      return 4;
    case ElementType::C3D10:
      return 10;
    case ElementType::C3D8:
    case ElementType::C3D8R:
      return 8;
    case ElementType::C3D20:
    case ElementType::C3D20R:
      return 20;
    case ElementType::C3D6:
      return 6;
    case ElementType::C3D15:
      return 15;
  }
  return 0;
}

std::string_view element_type_name(ElementType t);
bool parse_element_type(std::string_view name, ElementType& out);

struct Element {
  Index id{};
  ElementType type{ElementType::C3D4};
  std::vector<Index> nodes;  // node ids in CalculiX/Abaqus element ordering
};

}  // namespace cxpp
