#pragma once
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "calculixpp/core/element.hpp"
#include "calculixpp/core/node.hpp"

namespace cxpp {

// A named *SURFACE. TYPE=ELEMENT surfaces are a list of (element id, face 1..4)
// pairs; TYPE=NODE surfaces are a list of node ids. Exactly one list is populated
// per surface, selected by `type` (spec: mesh-and-model — *SURFACE).
struct Surface {
  enum class Type { Element, Node };
  std::string name;
  Type type{Type::Element};
  // (element id, face number 1..4) — populated when type == Element.
  std::vector<std::pair<Index, int>> faces;
  // node ids — populated when type == Node.
  std::vector<Index> nodes;
};

// In-memory finite element mesh: nodes, elements, and named sets/surfaces.
// (spec: mesh-and-model — Phase-1 subset: *NODE, *ELEMENT C3D4/C3D10,
// *NSET/*ELSET, *SURFACE.)
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

  // Store a named *SURFACE. A repeated name replaces the prior definition.
  void add_surface(Surface surface);
  const Surface* surface(const std::string& name) const;
  const std::unordered_map<std::string, Surface>& surfaces() const { return surfaces_; }

 private:
  std::vector<Node> nodes_;
  std::vector<Element> elements_;
  std::unordered_map<Index, Index> node_id_to_index_;
  std::unordered_map<Index, Index> element_id_to_index_;
  std::unordered_map<std::string, std::vector<Index>> nsets_;
  std::unordered_map<std::string, std::vector<Index>> elsets_;
  std::unordered_map<std::string, Surface> surfaces_;
};

}  // namespace cxpp
