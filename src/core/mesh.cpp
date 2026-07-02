#include "calculixpp/core/mesh.hpp"

#include <stdexcept>
#include <utility>

namespace cxpp {

Index Mesh::add_node(Index id, Vec3 x) {
  const auto idx = static_cast<Index>(nodes_.size());
  if (!node_id_to_index_.emplace(id, idx).second) {
    throw std::runtime_error("duplicate node id");
  }
  nodes_.push_back(Node{id, x});
  return idx;
}

Index Mesh::add_element(Index id, ElementType type, std::vector<Index> node_ids) {
  if (static_cast<int>(node_ids.size()) != nodes_per_element(type)) {
    throw std::runtime_error("element connectivity size mismatch");
  }
  const auto idx = static_cast<Index>(elements_.size());
  if (!element_id_to_index_.emplace(id, idx).second) {
    throw std::runtime_error("duplicate element id");
  }
  elements_.push_back(Element{id, type, std::move(node_ids)});
  return idx;
}

Index Mesh::node_index(Index id) const {
  const auto it = node_id_to_index_.find(id);
  return it == node_id_to_index_.end() ? Index{-1} : it->second;
}

Index Mesh::element_index(Index id) const {
  const auto it = element_id_to_index_.find(id);
  return it == element_id_to_index_.end() ? Index{-1} : it->second;
}

void Mesh::add_nset(const std::string& name, std::vector<Index> ids) {
  auto& set = nsets_[name];
  set.insert(set.end(), ids.begin(), ids.end());
}

void Mesh::add_elset(const std::string& name, std::vector<Index> ids) {
  auto& set = elsets_[name];
  set.insert(set.end(), ids.begin(), ids.end());
}

const std::vector<Index>* Mesh::nset(const std::string& name) const {
  const auto it = nsets_.find(name);
  return it == nsets_.end() ? nullptr : &it->second;
}

const std::vector<Index>* Mesh::elset(const std::string& name) const {
  const auto it = elsets_.find(name);
  return it == elsets_.end() ? nullptr : &it->second;
}

}  // namespace cxpp
