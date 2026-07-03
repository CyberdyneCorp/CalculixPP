#include "calculixpp/core/element.hpp"

#include <utility>

namespace cxpp {

std::string_view element_type_name(ElementType t) {
  switch (t) {
    case ElementType::C3D4:
      return "C3D4";
    case ElementType::C3D10:
      return "C3D10";
    case ElementType::C3D8:
      return "C3D8";
    case ElementType::C3D8R:
      return "C3D8R";
    case ElementType::C3D20:
      return "C3D20";
    case ElementType::C3D20R:
      return "C3D20R";
    case ElementType::C3D6:
      return "C3D6";
    case ElementType::C3D15:
      return "C3D15";
  }
  return "UNKNOWN";
}

bool parse_element_type(std::string_view name, ElementType& out) {
  static const std::pair<std::string_view, ElementType> kTable[] = {
      {"C3D4", ElementType::C3D4},     {"C3D10", ElementType::C3D10},
      {"C3D8", ElementType::C3D8},     {"C3D8R", ElementType::C3D8R},
      {"C3D20", ElementType::C3D20},   {"C3D20R", ElementType::C3D20R},
      {"C3D6", ElementType::C3D6},     {"C3D15", ElementType::C3D15},
  };
  for (const auto& [n, t] : kTable) {
    if (name == n) {
      out = t;
      return true;
    }
  }
  return false;
}

}  // namespace cxpp
