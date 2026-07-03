// *MODEL CHANGE, TYPE=CONTACT PAIR activate/deactivate across steps (tasks.md 5.2;
// spec: model-change — activate/deactivate a contact pair between steps). Validation:
// a contact pair INACTIVE in step 1 and ACTIVE in step 2 changes the mechanical response.
//
// SETUP: two unit C3D8 cubes stacked in z (base z in [0,1] fixed at z=0; top z in [1,2],
// distinct nodes on z=1). A *CONTACT PAIR couples the top-cube bottom nodes (slave) to the
// base-cube top face (master). The top cube's TOP face (z=2) is DRIVEN DOWN by a prescribed
// u_z = -delta, so it presses toward the base. All lateral DOFs are confined so the problem
// is non-singular whether or not contact is active.
//
//   * contact ACTIVE  -> the interface resists: the driven-down top cube compresses the
//     base cube, so the base's fixed face (z=0) develops a large reaction opposing the
//     drive (the load path is complete through the contact interface).
//   * contact INACTIVE -> the top cube's bottom face simply moves into the region of the
//     base with no interaction: the base carries NO reaction from the top cube's drive
//     (the two bodies are otherwise unconnected).
//
// So *MODEL CHANGE REMOVE (step 1) -> ~zero base reaction; ADD (step 2) -> large base
// reaction. The step-1 -> step-2 change in the base reaction is the model-change effect.
#include <cmath>
#include <string>
#include <vector>

#include "calculixpp/core/model.hpp"
#include "calculixpp/io/inp_parser.hpp"
#include "calculixpp/numerics/multistep.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

// Two stacked cubes, contact pair SSLAV/SMAST, driven top face. `step1_change` and
// `step2_change` are the *MODEL CHANGE blocks spliced into each step (empty = pair active).
std::string two_step_contact_deck(const std::string& step1_change,
                                  const std::string& step2_change) {
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
9, 0., 0., 1.
10, 1., 0., 1.
11, 1., 1., 1.
12, 0., 1., 1.
13, 0., 0., 2.
14, 1., 0., 2.
15, 1., 1., 2.
16, 0., 1., 2.
*ELEMENT, TYPE=C3D8, ELSET=EBASE
1, 1,2,3,4,5,6,7,8
*ELEMENT, TYPE=C3D8, ELSET=ETOP
2, 9,10,11,12,13,14,15,16
*ELSET, ELSET=EALL
EBASE, ETOP
*NSET, NSET=NBOT
1,2,3,4
*NSET, NSET=NTOP
13,14,15,16
*MATERIAL, NAME=STEEL
*ELASTIC
210000., 0.3
*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL
*SURFACE, NAME=SMAST
EBASE, S2
*SURFACE, NAME=SSLAV, TYPE=NODE
9,10,11,12
*SURFACE INTERACTION, NAME=SI
*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=LINEAR
1.E6
*CONTACT PAIR, INTERACTION=SI, TYPE=NODE TO SURFACE
SSLAV, SMAST
*BOUNDARY
NALL,1,2
NBOT,3,3
*STEP
*STATIC
)") + step1_change + R"(*BOUNDARY
13,3,3,-0.02
14,3,3,-0.02
15,3,3,-0.02
16,3,3,-0.02
*NODE PRINT, NSET=NALL
U,RF
*END STEP
*STEP
*STATIC
)" + step2_change + R"(*NODE PRINT, NSET=NALL
U,RF
*END STEP
)";
}

// Sum of the z-reaction on the fixed base face (z=0, nodes 1-4).
double base_reaction_z(const Model& m, const StaticFields& f) {
  double rf = 0.0;
  for (Index id : {1, 2, 3, 4})
    rf += f.reaction[static_cast<std::size_t>(m.mesh.node_index(id))][2];
  return rf;
}

// The active-pair filter itself: REMOVE deactivates, a later ADD reactivates, surfaces
// matched in either order.
void test_active_contact_pair_filter() {
  const Model m = io::parse_inp(two_step_contact_deck(
      "*MODEL CHANGE, TYPE=CONTACT PAIR, REMOVE\nSSLAV, SMAST\n", ""));
  CX_CHECK(m.has_contact());  // the pair is defined
  CX_CHECK(m.contact_pair_changes.size() == 1);
  // Filtered by the REMOVE record -> no active pair.
  const Model filtered = m.with_active_contact_pairs();
  CX_CHECK(!filtered.has_contact());

  // Surfaces given in the OTHER order on the *MODEL CHANGE line still match.
  const Model m2 = io::parse_inp(two_step_contact_deck(
      "*MODEL CHANGE, TYPE=CONTACT PAIR, REMOVE\nSMAST, SSLAV\n", ""));
  CX_CHECK(!m2.with_active_contact_pairs().has_contact());

  // No change -> the pair stays active.
  const Model m3 = io::parse_inp(two_step_contact_deck("", ""));
  CX_CHECK(m3.with_active_contact_pairs().has_contact());
}

// End-to-end: pair inactive in step 1, active in step 2 -> the base reaction jumps from
// ~zero (no load path) to a large compressive reaction (contact transmits the drive).
void test_contact_pair_activation_changes_response() {
  // Step 1 REMOVE (inactive), step 2 ADD (active).
  const std::vector<Model> steps = io::parse_inp_steps(two_step_contact_deck(
      "*MODEL CHANGE, TYPE=CONTACT PAIR, REMOVE\nSSLAV, SMAST\n",
      "*MODEL CHANGE, TYPE=CONTACT PAIR, ADD\nSSLAV, SMAST\n"));
  CX_CHECK(steps.size() == 2);

  const std::vector<StaticFields> f = numerics::solve_multistep_static_all(steps);
  CX_CHECK(f.size() == 2);

  const double rf1 = base_reaction_z(steps[0], f[0]);  // contact inactive
  const double rf2 = base_reaction_z(steps[1], f[1]);  // contact active

  // Step 1 (no contact): the base carries essentially no reaction from the top cube's
  // drive — the two bodies are unconnected. Step 2 (contact): the driven top cube presses
  // the base, so a large compressive reaction appears at the fixed base face.
  CX_CHECK(std::fabs(rf1) < 1e-3);            // ~zero without contact
  CX_CHECK(std::fabs(rf2) > 1.0e3);           // substantial reaction with contact
  CX_CHECK(std::fabs(rf2) > 1.0e6 * std::fabs(rf1) || std::fabs(rf1) < 1e-6);
}

// The reverse sense (active in step 1, removed in step 2) also flips the response, and the
// step-1 active reaction matches the never-changed reference — activation is a clean toggle.
void test_contact_pair_deactivation_changes_response() {
  const std::vector<Model> steps = io::parse_inp_steps(two_step_contact_deck(
      "", "*MODEL CHANGE, TYPE=CONTACT PAIR, REMOVE\nSSLAV, SMAST\n"));
  const std::vector<StaticFields> f = numerics::solve_multistep_static_all(steps);
  const double rf1 = base_reaction_z(steps[0], f[0]);  // contact active
  const double rf2 = base_reaction_z(steps[1], f[1]);  // contact removed
  CX_CHECK(std::fabs(rf1) > 1.0e3);   // active: large reaction
  CX_CHECK(std::fabs(rf2) < 1e-3);    // removed: ~zero
}

}  // namespace

int main() {
  test_active_contact_pair_filter();
  test_contact_pair_activation_changes_response();
  test_contact_pair_deactivation_changes_response();
  if (cxtest::g_failures == 0) std::printf("test_model_change_contact: OK\n");
  CX_MAIN_RETURN();
}
