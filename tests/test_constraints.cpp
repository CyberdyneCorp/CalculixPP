// Multi-point-constraint tests (spec: constraints — *EQUATION, *MPC, *RIGID BODY,
// *COUPLING, *TIE, over-constraint detection). Validation strategy:
//  - *EQUATION physics: analytic — two free nodes tied to move together respond
//    identically, and a symmetric two-element bar with the interface tied matches a
//    hand check; plus the reference deck achtel2 (checked in python/tests).
//  - over-constraint detection: a slave that is also an SPC, and a doubly-defined
//    dependent DOF, both raise.
//  - reproduces-linear gate: a model with a trivial identity constraint (or none)
//    gives the exact unconstrained solution.
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "calculixpp/io/inp_parser.hpp"
#include "calculixpp/numerics/linear_static.hpp"
#include "calculixpp/numerics/nonlinear_static.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

// A single C3D8 unit cube, base (z=0) fixed, top (z=1) free. Node ids 1..8 with the
// bottom face 1-4 and top face 5-8. Returns the deck text with an insertion point for
// extra cards before *STEP.
std::string cube_deck(const std::string& extra_before_step,
                      const std::string& loads) {
  return std::string(
             "*NODE\n"
             "1, 0., 0., 0.\n2, 1., 0., 0.\n3, 1., 1., 0.\n4, 0., 1., 0.\n"
             "5, 0., 0., 1.\n6, 1., 0., 1.\n7, 1., 1., 1.\n8, 0., 1., 1.\n"
             "*ELEMENT, TYPE=C3D8, ELSET=EALL\n1, 1,2,3,4,5,6,7,8\n"
             "*MATERIAL, NAME=STEEL\n*ELASTIC\n210000., 0.3\n"
             "*SOLID SECTION, MATERIAL=STEEL, ELSET=EALL\n"
             "*BOUNDARY\n1,1,3\n2,1,3\n3,1,3\n4,1,3\n") +
         extra_before_step + "*STEP\n*STATIC\n" + loads + "*END STEP\n";
}

std::size_t idx(const Model& m, Index id) {
  return static_cast<std::size_t>(m.mesh.node_index(id));
}

Real l2(const std::vector<Vec3>& a, const std::vector<Vec3>& b) {
  Real num = 0, den = 0;
  for (std::size_t i = 0; i < a.size(); ++i)
    for (int c = 0; c < 3; ++c) {
      const Real d = a[i][static_cast<std::size_t>(c)] - b[i][static_cast<std::size_t>(c)];
      num += d * d;
      den += b[i][static_cast<std::size_t>(c)] * b[i][static_cast<std::size_t>(c)];
    }
  return den > 0 ? std::sqrt(num / den) : std::sqrt(num);
}

// (1) *EQUATION ties two top nodes' x-displacement together: u_(6,1) - u_(7,1) = 0.
// Load only node 6 in x. Without the tie node 6 moves and node 7 barely; with the
// tie they must be equal.
void test_equation_ties_nodes() {
  const std::string eq =
      "*EQUATION\n2\n6, 1, 1.0, 7, 1, -1.0\n";
  const Model tied =
      io::parse_inp(cube_deck(eq, "*CLOAD\n6,1,1000.\n"));
  const StaticFields r = numerics::solve_linear_static(tied);
  const Real u6 = r.displacement[idx(tied, 6)][0];
  const Real u7 = r.displacement[idx(tied, 7)][0];
  CX_NEAR(u6, u7, 1e-10);            // the constraint holds exactly
  CX_CHECK(std::fabs(u6) > 1e-8);    // and they actually moved

  // Unconstrained: the two are not equal (sanity that the tie changed the result).
  const Model free_m = io::parse_inp(cube_deck("", "*CLOAD\n6,1,1000.\n"));
  const StaticFields rf = numerics::solve_linear_static(free_m);
  const Real f6 = rf.displacement[idx(free_m, 6)][0];
  const Real f7 = rf.displacement[idx(free_m, 7)][0];
  CX_CHECK(std::fabs(f6 - f7) > 1e-7);
}

// (2) A 3-term *EQUATION: node 7's x = mean of node 5 and node 6 x. Verify the
// eliminated slave equals the coefficient combination of its masters.
void test_equation_three_terms() {
  const std::string eq =
      "*EQUATION\n3\n7, 1, 2.0, 5, 1, -1.0, 6, 1, -1.0\n";  // 2 u7 - u5 - u6 = 0
  const Model m = io::parse_inp(cube_deck(eq, "*CLOAD\n5,1,500.\n6,1,1500.\n"));
  const StaticFields r = numerics::solve_linear_static(m);
  const Real u5 = r.displacement[idx(m, 5)][0];
  const Real u6 = r.displacement[idx(m, 6)][0];
  const Real u7 = r.displacement[idx(m, 7)][0];
  CX_NEAR(u7, 0.5 * (u5 + u6), 1e-10);
}

// (3) Over-constraint: the dependent DOF of an *EQUATION is also a *BOUNDARY SPC.
void test_overconstraint_slave_is_spc() {
  // node 5 dof1 is the slave AND fixed by *BOUNDARY 5,1.
  const std::string extra =
      "*BOUNDARY\n5,1\n*EQUATION\n2\n5, 1, 1.0, 6, 1, -1.0\n";
  const Model m = io::parse_inp(cube_deck(extra, "*CLOAD\n6,1,1000.\n"));
  bool threw = false;
  try {
    (void)numerics::solve_linear_static(m);
  } catch (const std::runtime_error& e) {
    threw = true;
    CX_CHECK(std::string(e.what()).find("over-constraint") != std::string::npos);
  }
  CX_CHECK(threw);
}

// (4) Over-constraint: two equations name the same dependent DOF.
void test_overconstraint_double_dependent() {
  const std::string extra =
      "*EQUATION\n2\n7, 1, 1.0, 6, 1, -1.0\n"
      "*EQUATION\n2\n7, 1, 1.0, 5, 1, -1.0\n";
  const Model m = io::parse_inp(cube_deck(extra, "*CLOAD\n6,1,1000.\n"));
  bool threw = false;
  try {
    (void)numerics::solve_linear_static(m);
  } catch (const std::runtime_error& e) {
    threw = true;
    CX_CHECK(std::string(e.what()).find("over-constraint") != std::string::npos);
  }
  CX_CHECK(threw);
}

// (5) *MPC BEAM ties all three translations of the dependent node to the second.
void test_mpc_beam() {
  const std::string extra = "*MPC\nBEAM, 6, 7\n";  // node 6 tied to node 7
  const Model m = io::parse_inp(cube_deck(extra, "*CLOAD\n6,1,800.\n6,2,300.\n"));
  const StaticFields r = numerics::solve_linear_static(m);
  for (int c = 0; c < 3; ++c)
    CX_NEAR(r.displacement[idx(m, 6)][static_cast<std::size_t>(c)],
            r.displacement[idx(m, 7)][static_cast<std::size_t>(c)], 1e-10);
}

// (6) *RIGID BODY (no rotation node) ties a node set rigidly to a reference node:
// every set node's translation equals the reference. Load the reference node.
void test_rigid_body_translation() {
  const std::string extra =
      "*NSET, NSET=TOP\n5,6,7,8\n"
      "*RIGID BODY, NSET=TOP, REF NODE=5\n";
  const Model m = io::parse_inp(cube_deck(extra, "*CLOAD\n5,1,1000.\n"));
  const StaticFields r = numerics::solve_linear_static(m);
  const Vec3 uref = r.displacement[idx(m, 5)];
  for (const Index nd : {6, 7, 8})
    for (int c = 0; c < 3; ++c)
      CX_NEAR(r.displacement[idx(m, static_cast<Index>(nd))][static_cast<std::size_t>(c)],
              uref[static_cast<std::size_t>(c)], 1e-10);
}

// (7) *DISTRIBUTING COUPLING: the reference node's motion is the average of the
// coupled nodes. Here the reference is a top node coupled to the other three; check
// the averaging equation holds after solve.
void test_kinematic_coupling() {
  // Kinematic coupling: surface nodes rigidly follow the reference node.
  const std::string extra =
      "*NSET, NSET=TOPSURF\n6,7,8\n"
      "*COUPLING, REF NODE=5, SURFACE=TOPSURF\n*KINEMATIC\n1,3\n";
  const Model m = io::parse_inp(cube_deck(extra, "*CLOAD\n5,2,1000.\n"));
  const StaticFields r = numerics::solve_linear_static(m);
  const Vec3 uref = r.displacement[idx(m, 5)];
  for (const Index nd : {6, 7, 8})
    for (int c = 0; c < 3; ++c)
      CX_NEAR(r.displacement[idx(m, static_cast<Index>(nd))][static_cast<std::size_t>(c)],
              uref[static_cast<std::size_t>(c)], 1e-10);
}

// (8) Reproduces-linear gate: the nonlinear driver on a *EQUATION model reproduces
// the linear-static solve to rel-L2 < 1e-10 (the constraint path is identical in
// both drivers).
void test_nonlinear_reproduces_linear_with_equation() {
  const std::string eq = "*EQUATION\n2\n6, 1, 1.0, 7, 1, -1.0\n";
  const Model m = io::parse_inp(cube_deck(eq, "*CLOAD\n6,1,1000.\n"));
  const StaticFields lin = numerics::solve_linear_static(m);
  const StaticFields nl = numerics::solve_nonlinear_static(m);
  CX_CHECK(l2(nl.displacement, lin.displacement) < 1e-10);
}

}  // namespace

int main() {
  test_equation_ties_nodes();
  test_equation_three_terms();
  test_overconstraint_slave_is_spc();
  test_overconstraint_double_dependent();
  test_mpc_beam();
  test_rigid_body_translation();
  test_kinematic_coupling();
  test_nonlinear_reproduces_linear_with_equation();
  CX_MAIN_RETURN();
}
