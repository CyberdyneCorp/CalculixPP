#pragma once
#include <string>
#include <vector>

#include "calculixpp/core/types.hpp"

namespace cxpp {

// *MPC: an analytic multi-point constraint expanded into linear equations.
// (spec: constraints — Multi-point constraints.)
//   PLANE:    the listed nodes are kept in the plane defined by the first three
//             (dependent) nodes' motion; each extra node contributes one equation.
//   STRAIGHT: the listed nodes are kept on the straight line through the first two.
//   BEAM:     a rigid beam between the first (dependent) node and the second node
//             (all three translational DOFs tied — small-displacement rigid link).
struct Mpc {
  enum class Kind { Plane, Straight, Beam, User };
  Kind kind{Kind::Beam};
  std::vector<Index> nodes;  // node ids in card order (first is dependent)
  std::string user_name;     // for Kind::User: the registered hook name
};

// *RIGID BODY: a node set tied to a reference (and optional rotation) node so the set
// moves as a rigid body. In small displacement each set node's translation equals the
// reference translation plus rot x (x_node - x_ref). Without a rotation node the set
// is tied purely to the reference translation. (spec: constraints — Rigid bodies.)
struct RigidBody {
  std::string nset;       // node set that moves rigidly (empty -> use `nodes`)
  std::vector<Index> nodes;
  Index ref_node{0};      // reference node (translation)
  Index rot_node{0};      // rotation node (0 -> none; rotational DOFs stored on it)
};

// *COUPLING / *DISTRIBUTING COUPLING: couples a surface's nodes to a reference node.
// KINEMATIC constrains the surface rigidly to the reference (each surface node's
// translation = ref translation, small-displacement); DISTRIBUTING spreads the
// reference node's motion/load across the surface by weights (average constraint).
// (spec: constraints — Distributing and kinematic couplings.)
struct Coupling {
  enum class Kind { Kinematic, Distributing };
  Kind kind{Kind::Kinematic};
  Index ref_node{0};
  std::string surface;             // surface / node-set name the coupling acts on
  std::vector<Index> nodes;        // resolved surface nodes (filled by the parser)
  std::vector<int> dofs{1, 2, 3};  // constrained DOFs (default translations)
};

// *TIE: bond two surfaces. For matching meshes each slave node is tied to the
// coincident master node (all translations equal). Non-matching (mortar) support is
// partial. (spec: constraints — Surface ties.)
struct Tie {
  std::string name;
  std::string slave_surface;
  std::string master_surface;
  std::vector<Index> slave_nodes;   // resolved (filled by the parser)
  std::vector<Index> master_nodes;  // resolved (filled by the parser)
  Real tolerance{0.0};              // position tolerance for node association
};

}  // namespace cxpp
