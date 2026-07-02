#pragma once
#include <string_view>
#include <vector>

#include "calculixpp/core/types.hpp"

namespace cxpp {

// Phase-1 element library: linear and quadratic tetrahedra.
enum class ElementType { C3D4, C3D10 };

constexpr int nodes_per_element(ElementType t) {
  return t == ElementType::C3D4 ? 4 : 10;
}

std::string_view element_type_name(ElementType t);
bool parse_element_type(std::string_view name, ElementType& out);

struct Element {
  Index id{};
  ElementType type{ElementType::C3D4};
  std::vector<Index> nodes;  // node ids in CalculiX/Abaqus element ordering
};

}  // namespace cxpp
