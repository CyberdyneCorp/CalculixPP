#pragma once
#include <string>

namespace cxpp {

// *SOLID SECTION binding a material to an element set (spec: element-sections — Phase 1).
struct SolidSection {
  std::string elset;
  std::string material;
};

}  // namespace cxpp
