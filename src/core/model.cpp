#include "calculixpp/core/model.hpp"

#include <stdexcept>
#include <string>

namespace cxpp {

std::vector<ElasticIso> Model::element_elastic() const {
  std::vector<ElasticIso> out(mesh.num_elements());
  std::vector<bool> assigned(mesh.num_elements(), false);

  for (const auto& section : sections) {
    const std::vector<Index>* ids = mesh.elset(section.elset);
    if (ids == nullptr) {
      throw std::runtime_error("solid section references unknown elset '" +
                               section.elset + "'");
    }
    const auto mat = materials.find(section.material);
    if (mat == materials.end()) {
      throw std::runtime_error("solid section references unknown material '" +
                               section.material + "'");
    }
    if (!mat->second.elastic) {
      throw std::runtime_error("material '" + section.material +
                               "' has no *ELASTIC data");
    }
    for (const Index elem_id : *ids) {
      const Index ei = mesh.element_index(elem_id);
      if (ei < 0) {
        throw std::runtime_error("elset '" + section.elset +
                                 "' references unknown element");
      }
      out[static_cast<std::size_t>(ei)] = *mat->second.elastic;
      assigned[static_cast<std::size_t>(ei)] = true;
    }
  }

  for (std::size_t i = 0; i < assigned.size(); ++i) {
    if (!assigned[i]) {
      throw std::runtime_error("element " +
                               std::to_string(mesh.elements()[i].id) +
                               " has no solid section / material");
    }
  }
  return out;
}

}  // namespace cxpp
