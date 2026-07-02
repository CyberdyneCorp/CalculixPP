#include "calculixpp/core/element.hpp"

namespace cxpp {

std::string_view element_type_name(ElementType t) {
  switch (t) {
    case ElementType::C3D4:
      return "C3D4";
    case ElementType::C3D10:
      return "C3D10";
  }
  return "UNKNOWN";
}

bool parse_element_type(std::string_view name, ElementType& out) {
  if (name == "C3D4") {
    out = ElementType::C3D4;
    return true;
  }
  if (name == "C3D10") {
    out = ElementType::C3D10;
    return true;
  }
  return false;
}

}  // namespace cxpp
