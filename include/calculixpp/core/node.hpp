#pragma once
#include "calculixpp/core/types.hpp"

namespace cxpp {

struct Node {
  Index id{};   // user-facing node id from the deck
  Vec3 x{};     // coordinates (missing coords default to 0, per Abaqus convention)
};

}  // namespace cxpp
