// End-to-end validation of the hex families (C3D8 full + C3D8R reduced) on a
// single-element uniaxial-tension patch with a known analytical solution.
//
// A unit cube [0,1]^3 is stretched in x by a prescribed displacement Delta on the
// x=1 face while the x=0 face is held at u_x=0. The problem is uniaxial: the exact
// solution is a homogeneous strain e_xx = Delta, e_yy = e_zz = -nu*Delta, all shear
// zero, with uniform stress s_xx = E*Delta and zero lateral stress. A single hex
// (full or reduced) must reproduce this exactly, and the x-face reaction must equal
// the analytical axial force s_xx * Area. This is the constant-stress patch test
// required for the new hex physics.
#include <array>
#include <cmath>
#include <vector>

#include "calculixpp/numerics/linear_static.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

// Unit cube in C3D8 node order (bottom face z=0 CCW, then top z=1).
const std::array<Vec3, 8> kCube{{
    {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
    {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1},
}};

// Nodes on the x=0 face (1,4,5,8) and x=1 face (2,3,6,7) — 1-based ids.
const std::array<Index, 4> kFaceX0{1, 4, 5, 8};
const std::array<Index, 4> kFaceX1{2, 3, 6, 7};

Model make_cube(ElementType type, Real E, Real nu, Real delta) {
  Model m;
  for (int k = 0; k < 8; ++k) m.mesh.add_node(k + 1, kCube[static_cast<std::size_t>(k)]);
  std::vector<Index> conn;
  for (int k = 0; k < 8; ++k) conn.push_back(k + 1);
  m.mesh.add_element(1, type, conn);
  m.mesh.add_elset("EALL", {1});
  m.materials["EL"] = Material{"EL", ElasticIso{E, nu}, std::nullopt};
  m.sections.push_back(SolidSection{"EALL", "EL"});

  // x=0 face: u_x = 0. Also pin y,z on one node (node 1) and z on node 4 to remove
  // the two lateral rigid-body translations and the three rotations without
  // over-constraining the Poisson contraction (statically determinate support).
  for (Index nd : kFaceX0) m.spcs.push_back(Spc{nd, 1, 0.0});
  m.spcs.push_back(Spc{1, 2, 0.0});  // node at origin: fix y
  m.spcs.push_back(Spc{1, 3, 0.0});  // and z
  m.spcs.push_back(Spc{4, 3, 0.0});  // node (0,1,0): fix z -> kills rot about x
  m.spcs.push_back(Spc{5, 2, 0.0});  // node (0,0,1): fix y -> kills rot about...
  // x=1 face: prescribe u_x = delta.
  for (Index nd : kFaceX1) m.spcs.push_back(Spc{nd, 1, delta});
  return m;
}

// Check the homogeneous uniaxial solution: e_xx = delta everywhere, lateral
// contraction -nu*delta, and the total x-reaction on the fixed face = E*delta*Area.
void run_uniaxial(ElementType type) {
  const Real E = 210000.0, nu = 0.3, delta = 0.01;
  const Model m = make_cube(type, E, nu, delta);
  const StaticFields res = numerics::solve_linear_static(m);

  // The x=1 face nodes moved by exactly delta in x (prescribed) and the free faces
  // contract by -nu*delta in y/z. Check a representative free node: node 7 (1,1,1).
  const Index i7 = m.mesh.node_index(7);
  CX_NEAR(res.displacement[static_cast<std::size_t>(i7)][0], delta, 1e-10);
  CX_NEAR(res.displacement[static_cast<std::size_t>(i7)][1], -nu * delta, 1e-8);
  CX_NEAR(res.displacement[static_cast<std::size_t>(i7)][2], -nu * delta, 1e-8);

  // Nodal stress must be uniform uniaxial: s_xx = E*delta, lateral ~ 0.
  const Real sxx = res.stress[static_cast<std::size_t>(i7)][0];
  CX_NEAR(sxx, E * delta, 1e-4 * E * delta);
  CX_NEAR(res.stress[static_cast<std::size_t>(i7)][1], 0.0, 1e-4 * E * delta);
  CX_NEAR(res.stress[static_cast<std::size_t>(i7)][2], 0.0, 1e-4 * E * delta);

  // Total reaction in x on the x=0 face equals the applied axial force = s_xx*Area
  // (Area = 1). Newton's third law: sum of x-reactions on x=1 face = -that.
  Real rf0 = 0.0;
  for (Index nd : kFaceX0)
    rf0 += res.reaction[static_cast<std::size_t>(m.mesh.node_index(nd))][0];
  CX_NEAR(rf0, -E * delta, 1e-4 * E * delta);
}

}  // namespace

int main() {
  run_uniaxial(ElementType::C3D8);
  run_uniaxial(ElementType::C3D8R);
  if (cxtest::g_failures == 0) std::printf("test_hex: OK\n");
  CX_MAIN_RETURN();
}
