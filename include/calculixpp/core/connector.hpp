#pragma once
#include <vector>

#include "calculixpp/core/types.hpp"

// Discrete / connector elements (spec: element-sections — *MASS, *SPRING,
// *DASHPOT). These are NOT isoparametric solids: they carry no shape functions or
// integration rule. A connector references one or two mesh nodes and a scalar
// property, and contributes directly to the assembled system (springs) or is stored
// for later dynamics (mass, dashpot).
namespace cxpp {

// Spring connector. Three CalculiX forms, distinguished by `kind`:
//   Grounded : SPRING1 — one node, one DOF, grounded (reaction against ground).
//   Dof      : SPRING2 — two nodes, a DOF on each; force ~ k*(u1_dof1 - u2_dof2).
//   Axial    : SPRINGA — two nodes; force acts along the line joining them, so the
//              3x3 nodal block is k * (n outer n) with n the unit direction.
struct Spring {
  enum class Kind { Grounded, Dof, Axial };
  Kind kind{Kind::Axial};
  Index node1{};       // first (or only) node id
  Index node2{};       // second node id (Dof/Axial); unused for Grounded
  int dof1{1};         // acting DOF on node1 (Grounded/Dof); 1..3
  int dof2{1};         // acting DOF on node2 (Dof); 1..3
  Real stiffness{0.0}; // linear spring constant k
};

// Point mass at a node (*MASS + TYPE=MASS element). Stored for dynamics (mass
// matrix); contributes nothing to the static stiffness.
struct PointMass {
  Index node{};
  Real mass{0.0};
};

// Viscous dashpot connector (*DASHPOT + TYPE=DASHPOTA/DASHPOT1/DASHPOT2). Stored
// for dynamics (damping matrix); contributes nothing to the static stiffness.
struct Dashpot {
  enum class Kind { Grounded, Dof, Axial };
  Kind kind{Kind::Axial};
  Index node1{};
  Index node2{};
  int dof1{1};
  int dof2{1};
  Real coefficient{0.0};  // viscous damping coefficient c
};

}  // namespace cxpp
