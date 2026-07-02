#pragma once
#include <string>
#include <unordered_map>
#include <vector>

#include "calculixpp/core/element.hpp"
#include "calculixpp/core/node.hpp"

namespace cxpp {

// In-memory finite element mesh: nodes, elements, and named sets.
// (spec: mesh-and-model — Phase-1 subset: *NODE, *ELEMENT C3D4/C3D10, *NSET/*ELSET.)
class Mesh {
 public:
  // Returns the internal (dense) index of the added entity.
  Index add_node(Index id, Vec3 x);
  Index add_element(Index id, ElementType type, std::vector<Index> node_ids);

  const std::vector<Node>& nodes() const { return nodes_; }
  const std::vector<Element>& elements() const { return elements_; }
  std::size_t num_nodes() const { return nodes_.size(); }
  std::size_t num_elements() const { return elements_.size(); }

  // id -> internal index, or -1 if absent.
  Index node_index(Index id) const;
  Index element_index(Index id) const;

  void add_nset(const std::string& name, std::vector<Index> ids);
  void add_elset(const std::string& name, std::vector<Index> ids);
  const std::vector<Index>* nset(const std::string& name) const;
  const std::vector<Index>* elset(const std::string& name) const;

 private:
  std::vector<Node> nodes_;
  std::vector<Element> elements_;
  std::unordered_map<Index, Index> node_id_to_index_;
  std::unordered_map<Index, Index> element_id_to_index_;
  std::unordered_map<std::string, std::vector<Index>> nsets_;
  std::unordered_map<std::string, std::vector<Index>> elsets_;
};

}  // namespace cxpp
