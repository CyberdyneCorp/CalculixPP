#pragma once
#include <string>
#include <vector>

#include "calculixpp/core/element.hpp"
#include "calculixpp/core/types.hpp"

namespace cxpp {

// Global-analysis result carrier for submodeling (spec: submodeling — global result
// source). A `GlobalSolution` is a lightweight in-memory view of the global (coarse)
// mesh and its nodal displacement field: enough to locate the host element of any
// point and interpolate the displacement there with the element shape functions. It
// is deliberately source-agnostic — built directly in tests/bindings here; reading it
// from a global `.frd`/results file on disk is a documented follow-on.
//
// `elem_type`/`elem_conn` describe the global elements: `elem_conn[e]` lists the
// 0-based INTERNAL node indices (into `node_ids`/`coords`/`displacement`) of element e,
// in CalculiX/Abaqus element ordering. `displacement[i]` is the global nodal
// displacement U at global node i (aligned with `node_ids`/`coords`).
struct GlobalSolution {
  std::vector<Index> node_ids;               // global node ids (index i is the node's slot)
  std::vector<Vec3> coords;                  // global node coordinates (aligned with node_ids)
  std::vector<Vec3> displacement;            // global nodal displacement U (aligned with node_ids)
  std::vector<ElementType> elem_type;        // per-element type
  std::vector<std::vector<Index>> elem_conn; // per-element internal node indices (element order)

  std::size_t num_nodes() const { return node_ids.size(); }
  std::size_t num_elements() const { return elem_type.size(); }
};

// A parsed `*SUBMODEL` card (spec: submodeling — boundary declaration). The
// displacement-driven `TYPE=NODE` slice: `boundary_nset` is the driven boundary node
// set; `global_elset` (optional) bounds the global element search — empty means every
// global element is searched. `TYPE=SURFACE` is deferred.
struct SubmodelSpec {
  std::string boundary_nset;  // driven boundary node set name
  std::string global_elset;   // global element set searched for host elements (empty -> all)
};

}  // namespace cxpp
