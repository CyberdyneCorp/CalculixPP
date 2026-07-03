// *MODEL CHANGE, TYPE=ELEMENT element birth-death (tasks.md 5.1) and
// TYPE=CONTACT PAIR parse-only storage (tasks.md 5.2).
//
// Physics validation is by an ANALYTICAL check on a two-element uniaxial bar:
// two unit cubes stacked in y (sharing the y=1 face), each independently
// supported on its x=0 face and stretched by a prescribed u_x = delta on its
// x=1 face. Every element is a pure uniaxial specimen, so the total x-reaction
// on the fixed face is E*delta*Area with Area = number of ACTIVE cubes (each
// unit-area). This makes element removal a clean, closed-form softening test:
//
//   * both active      -> RF_x = 2 E delta   (Area = 2)
//   * one removed      -> RF_x = 1 E delta   (Area = 1; the removed cube carries
//                                             no stiffness, stress, or reaction)
//   * removed then ADDed back -> RF_x = 2 E delta again, IDENTICALLY equal to the
//                                never-removed model (strain-free reactivation:
//                                the single-step solve starts undeformed).
#include <cmath>
#include <string>

#include "calculixpp/io/inp_parser.hpp"
#include "calculixpp/numerics/linear_static.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

constexpr Real kE = 210000.0;
constexpr Real kNu = 0.3;
constexpr Real kDelta = 0.01;

// Confined (oedometer) modulus M = E(1-nu)/((1+nu)(1-2nu)): with all lateral DOFs
// held at zero, the axial stress is s_xx = M * e_xx, so the x-reaction resisting a
// prescribed u_x = delta on a unit-area cube is -M * delta (negative: it opposes
// the +x pull). A stack of N active unit cubes gives -N * M * delta.
constexpr Real kM = kE * (1.0 - kNu) / ((1.0 + kNu) * (1.0 - 2.0 * kNu));
// Per-active-cube fixed-face reaction (unit area).
constexpr Real kRfPerCube = -kM * kDelta;

// Confined uniaxial: every node has u_y = u_z = 0 (lateral confinement), the x=0
// face (nodes 1,4,5,8,9,11) has u_x = 0, and the x=1 face (2,3,6,7,10,12) has
// u_x = delta. Confining every node means NO node has a free DOF, so removing an
// element never orphans a DOF (the removed cube's exclusive nodes 9-12 stay fully
// prescribed) — the removal test stays non-singular while still being a clean
// analytical specimen: each unit cube carries the confined-modulus axial reaction.
std::string x0_and_pins() {
  std::string s;
  for (int nd = 1; nd <= 12; ++nd) {
    s += std::to_string(nd) + ", 2, 2\n";   // u_y = 0
    s += std::to_string(nd) + ", 3, 3\n";   // u_z = 0
  }
  for (Index nd : {1, 4, 5, 8, 9, 11})
    s += std::to_string(nd) + ", 1, 1\n";    // x=0 face: u_x = 0
  return s;
}

// x=1 faces of both cubes prescribed u_x = delta (nodes 2,3,6,7,10,12).
std::string prescribed_x1() {
  return R"(2, 1, 1, 0.01
3, 1, 1, 0.01
6, 1, 1, 0.01
7, 1, 1, 0.01
10, 1, 1, 0.01
12, 1, 1, 0.01
)";
}

// Two unit C3D8 cubes stacked in y: cube 1 spans y in [0,1], cube 2 y in [1,2],
// sharing the y=1 face (nodes 4,3,8,7). `model_change` is spliced in before
// *END STEP (empty for the baseline). Each cube's x=0 face is fixed (u_x=0) and
// its x=1 face prescribed u_x=delta, plus a statically-determinate set of pins to
// remove rigid-body motion without fighting the Poisson contraction.
std::string bar_deck(const std::string& model_change) {
  return std::string(R"(
*NODE, NSET=NALL
1, 0., 0., 0.
2, 1., 0., 0.
3, 1., 1., 0.
4, 0., 1., 0.
5, 0., 0., 1.
6, 1., 0., 1.
7, 1., 1., 1.
8, 0., 1., 1.
9,  0., 2., 0.
10, 1., 2., 0.
11, 0., 2., 1.
12, 1., 2., 1.
*ELEMENT, TYPE=C3D8, ELSET=CUBE1
1, 1, 2, 3, 4, 5, 6, 7, 8
*ELEMENT, TYPE=C3D8, ELSET=CUBE2
2, 4, 3, 10, 9, 8, 7, 12, 11
*ELSET, ELSET=EALL
1, 2
*MATERIAL, NAME=EL
*ELASTIC
210000., 0.3
*SOLID SECTION, ELSET=EALL, MATERIAL=EL
*BOUNDARY
)") + x0_and_pins() + prescribed_x1() +
         "*STEP\n*STATIC\n" + model_change + "*END STEP\n";
}

// Total x-reaction on the x=0 face of the model (sum over the fixed-face nodes).
Real fixed_face_reaction_x(const Model& m, const StaticFields& res) {
  Real rf = 0.0;
  for (Index nd : {1, 4, 5, 8, 9, 11})
    rf += res.reaction[static_cast<std::size_t>(m.mesh.node_index(nd))][0];
  return rf;
}

// Baseline: both cubes active -> RF_x = 2 * (-M delta) (two unit-area specimens).
void test_full_bar_reaction() {
  const Model m = io::parse_inp(bar_deck(""));
  CX_CHECK(m.deactivated_elements.empty());
  const StaticFields res = numerics::solve_linear_static(m);
  CX_NEAR(fixed_face_reaction_x(m, res), 2.0 * kRfPerCube, 1e-4 * std::fabs(kRfPerCube));
}

// REMOVE cube 2 -> it carries no stiffness/stress/reaction, so RF_x halves to one
// cube's worth. This is the "removing an element softens / creates a hole" check.
void test_removed_element_carries_no_load() {
  const Model m = io::parse_inp(
      bar_deck("*MODEL CHANGE, TYPE=ELEMENT, REMOVE\n2\n"));
  CX_CHECK(m.deactivated_elements.size() == 1);
  CX_CHECK(m.deactivated_elements[0] == 2);
  const auto mask = m.element_active_mask();
  CX_CHECK(mask[static_cast<std::size_t>(m.mesh.element_index(1))]);
  CX_CHECK(!mask[static_cast<std::size_t>(m.mesh.element_index(2))]);

  const StaticFields res = numerics::solve_linear_static(m);
  CX_NEAR(fixed_face_reaction_x(m, res), 1.0 * kRfPerCube, 1e-4 * std::fabs(kRfPerCube));

  // Cube-2-only interior node 10 sits on cube 2 alone; with cube 2 removed it
  // receives no stress from any active element.
  const Index i10 = m.mesh.node_index(10);
  CX_NEAR(res.stress[static_cast<std::size_t>(i10)][0], 0.0, 1e-6 * std::fabs(kRfPerCube));
}

// REMOVE then ADD the same set -> back to full model, IDENTICALLY. Reactivation
// is strain-free: the single-step solve starts undeformed, so the re-added cube
// contributes exactly as in the never-removed model (bit-identical reaction).
void test_reactivation_is_strain_free() {
  const Model full = io::parse_inp(bar_deck(""));
  const Model readded = io::parse_inp(bar_deck(
      "*MODEL CHANGE, TYPE=ELEMENT, REMOVE\n2\n"
      "*MODEL CHANGE, TYPE=ELEMENT, ADD\n2\n"));
  CX_CHECK(readded.deactivated_elements.empty());  // ADD cleared the REMOVE

  const StaticFields rf = numerics::solve_linear_static(full);
  const StaticFields rr = numerics::solve_linear_static(readded);
  // Every nodal displacement and reaction matches to solver precision.
  for (std::size_t i = 0; i < full.mesh.num_nodes(); ++i)
    for (int c = 0; c < 3; ++c) {
      CX_NEAR(rr.displacement[i][static_cast<std::size_t>(c)],
              rf.displacement[i][static_cast<std::size_t>(c)], 1e-12);
      CX_NEAR(rr.reaction[i][static_cast<std::size_t>(c)],
              rf.reaction[i][static_cast<std::size_t>(c)], 1e-9);
    }
}

// An elset name (not just a bare id) resolves on the *MODEL CHANGE data line.
void test_remove_by_elset_name() {
  const Model m = io::parse_inp(
      bar_deck("*MODEL CHANGE, TYPE=ELEMENT, REMOVE\nCUBE2\n"));
  CX_CHECK(m.deactivated_elements.size() == 1);
  CX_CHECK(m.deactivated_elements[0] == 2);
}

// A deck with no *MODEL CHANGE leaves every element active (byte-for-byte the
// pre-model-change mechanical path).
void test_no_model_change_all_active() {
  const Model m = io::parse_inp(bar_deck(""));
  const auto mask = m.element_active_mask();
  for (bool a : mask) CX_CHECK(a);
}

// TYPE=CONTACT PAIR is parsed and stored (consumed by the contact workflow),
// carrying the two surface names and the ADD/REMOVE sense.
void test_contact_pair_stored() {
  const Model m = io::parse_inp(bar_deck(
      "*MODEL CHANGE, TYPE=CONTACT PAIR, REMOVE\nSURF_A, SURF_B\n"));
  CX_CHECK(m.contact_pair_changes.size() == 1);
  CX_CHECK(m.contact_pair_changes[0].surface_a == "SURF_A");
  CX_CHECK(m.contact_pair_changes[0].surface_b == "SURF_B");
  CX_CHECK(!m.contact_pair_changes[0].add);  // REMOVE
  CX_CHECK(m.deactivated_elements.empty());  // contact pair does not touch elements
}

// Neither (or both) of ADD/REMOVE is an error.
void test_add_remove_required() {
  bool threw = false;
  try {
    io::parse_inp(bar_deck("*MODEL CHANGE, TYPE=ELEMENT\n2\n"));
  } catch (const io::ParseError&) {
    threw = true;
  }
  CX_CHECK(threw);
}

// An unknown TYPE is rejected.
void test_unknown_type_throws() {
  bool threw = false;
  try {
    io::parse_inp(bar_deck("*MODEL CHANGE, TYPE=NODE, REMOVE\n2\n"));
  } catch (const io::ParseError&) {
    threw = true;
  }
  CX_CHECK(threw);
}

}  // namespace

int main() {
  test_full_bar_reaction();
  test_removed_element_carries_no_load();
  test_reactivation_is_strain_free();
  test_remove_by_elset_name();
  test_no_model_change_all_active();
  test_contact_pair_stored();
  test_add_remove_required();
  test_unknown_type_throws();
  if (cxtest::g_failures == 0) std::printf("test_model_change: OK\n");
  CX_MAIN_RETURN();
}
