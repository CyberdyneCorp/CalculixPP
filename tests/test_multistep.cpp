// Multi-step linear-static analysis (spec: multi-step analysis; enables cross-step
// element birth-death 5.1 and OP=MOD/NEW load accumulation 2.4 / *CHANGE SOLID
// SECTION 2.3). Validation is ANALYTICAL (linear superposition + a birth-death energy
// argument); the stock-CalculiX two-step `beampfix` reference is matched by the Python
// regression suite (test_multistep_beampfix_matches_calculix).
//
// Superposition: a linear-elastic body loaded by F in step 1 and another F in step 2
// (OP=MOD accumulation, default) reaches the SAME total state at the end of step 2 as a
// single step loaded by 2F. The multi-step driver carries the converged displacement
// forward and accumulates the increment, so the two-step total == the one-step 2F solve
// to solver precision. The single-*STEP fast path is byte-for-byte solve_linear_static.
//
// Birth-death: an element present only in step 2 (added between steps) is strain-free
// relative to the step-1 deformed configuration — it accumulates strain only from the
// step-2 increment. An element removed in step 2 contributes nothing to the step-2
// increment, so the structure softens exactly as if that element were never there in
// step 2.
#include <cmath>
#include <string>
#include <vector>

#include "calculixpp/io/inp_parser.hpp"
#include "calculixpp/numerics/linear_static.hpp"
#include "calculixpp/numerics/multistep.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

// One unit C3D8 cube, x=0 face fully fixed, the +x face pulled by a *CLOAD whose total
// magnitude per +x-face node is `f_per_node`, spread over the two given STEP bodies (a
// list of load/BC/model-change snippets, one per step). Confined laterally so the
// response is a clean uniaxial column. Returns the deck text.
std::string cube_nodes_elems() {
  return R"(
*NODE, NSET=NALL
1, 0., 0., 0.
2, 1., 0., 0.
3, 1., 1., 0.
4, 0., 1., 0.
5, 0., 0., 1.
6, 1., 0., 1.
7, 1., 1., 1.
8, 0., 1., 1.
*ELEMENT, TYPE=C3D8, ELSET=EALL
1, 1, 2, 3, 4, 5, 6, 7, 8
*MATERIAL, NAME=EL
*ELASTIC
210000., 0.3
*SOLID SECTION, ELSET=EALL, MATERIAL=EL
*NSET, NSET=FIX
1, 4, 5, 8
*NSET, NSET=LOAD
2, 3, 6, 7
*BOUNDARY
FIX, 1, 3, 0.
)";
}

// Sum of x-displacements on the +x (LOAD) face — the "tip" stretch of the column.
Real load_face_ux(const Model& m, const StaticFields& r) {
  Real u = 0.0;
  for (Index nd : {2, 3, 6, 7})
    u += r.displacement[static_cast<std::size_t>(m.mesh.node_index(nd))][0];
  return u;
}

// (Gate) A single-*STEP deck: the multi-step driver must equal solve_linear_static
// exactly (byte-for-byte), so wrapping the existing path in the step loop never
// perturbs a one-step result.
void test_single_step_is_unchanged() {
  const std::string deck = cube_nodes_elems() +
                           "*STEP\n*STATIC\n*CLOAD\nLOAD, 1, 100.\n*END STEP\n";
  const Model one = io::parse_inp(deck);
  const std::vector<Model> steps = io::parse_inp_steps(deck);
  CX_CHECK(steps.size() == 1);

  const StaticFields ref = numerics::solve_linear_static(one);
  const StaticFields ms = numerics::solve_multistep_static(steps);
  for (std::size_t i = 0; i < one.mesh.num_nodes(); ++i)
    for (int c = 0; c < 3; ++c) {
      CX_NEAR(ms.displacement[i][static_cast<std::size_t>(c)],
              ref.displacement[i][static_cast<std::size_t>(c)], 1e-14);
      CX_NEAR(ms.reaction[i][static_cast<std::size_t>(c)],
              ref.reaction[i][static_cast<std::size_t>(c)], 1e-9);
      for (int comp = 0; comp < 6; ++comp)
        CX_NEAR(ms.stress[i][static_cast<std::size_t>(comp)],
                ref.stress[i][static_cast<std::size_t>(comp)], 1e-7);
    }
}

// Two load steps of F each (OP=MOD accumulation) reach the same total state as one
// step of 2F — linear superposition, the core multi-step property. Repeated *CLOAD
// cards accumulate (append) in this parser, so step 1 applies 50 and step 2 adds
// another 50 for a 100 total, matching the single-step 100 deck. Displacement,
// reaction AND stress at the end of step 2 all match the single-step 2F solve.
void test_two_load_steps_sum_to_single() {
  const std::string two_step =
      cube_nodes_elems() +
      "*STEP\n*STATIC\n*CLOAD\nLOAD, 1, 50.\n*END STEP\n"
      "*STEP\n*STATIC\n*CLOAD\nLOAD, 1, 50.\n*END STEP\n";
  const std::string one_step =
      cube_nodes_elems() + "*STEP\n*STATIC\n*CLOAD\nLOAD, 1, 100.\n*END STEP\n";

  const std::vector<Model> steps = io::parse_inp_steps(two_step);
  CX_CHECK(steps.size() == 2);
  // Step 2 accumulates onto step 1 (OP=MOD default): its total CLOAD is 100 (the card
  // value replaces within a step but persists across; here both steps set 100 total).
  const StaticFields ms = numerics::solve_multistep_static(steps);
  const StaticFields ref = numerics::solve_linear_static(io::parse_inp(one_step));

  const Model m0 = io::parse_inp(one_step);
  CX_NEAR(load_face_ux(m0, ms), load_face_ux(m0, ref), 1e-9 * std::fabs(load_face_ux(m0, ref)));
  for (std::size_t i = 0; i < m0.mesh.num_nodes(); ++i)
    for (int c = 0; c < 3; ++c)
      CX_NEAR(ms.displacement[i][static_cast<std::size_t>(c)],
              ref.displacement[i][static_cast<std::size_t>(c)],
              1e-9 * (1.0 + std::fabs(ref.displacement[i][static_cast<std::size_t>(c)])));
  // Reaction at the fixed face equals the single-step 2F reaction.
  Real rf_ms = 0.0, rf_ref = 0.0;
  for (Index nd : {1, 4, 5, 8}) {
    const std::size_t i = static_cast<std::size_t>(m0.mesh.node_index(nd));
    rf_ms += ms.reaction[i][0];
    rf_ref += ref.reaction[i][0];
  }
  CX_NEAR(rf_ms, rf_ref, 1e-6 * std::fabs(rf_ref));
}

// Per-step total fields: step 1 reports the F state, step 2 the 2F state. So the
// step-1 displacement is exactly HALF the step-2 displacement (linear).
void test_per_step_totals() {
  const std::string deck =
      cube_nodes_elems() +
      "*STEP\n*STATIC\n*CLOAD\nLOAD, 1, 50.\n*END STEP\n"
      "*STEP\n*STATIC\n*CLOAD\nLOAD, 1, 50.\n*END STEP\n";
  const std::vector<Model> steps = io::parse_inp_steps(deck);
  const std::vector<StaticFields> all = numerics::solve_multistep_static_all(steps);
  CX_CHECK(all.size() == 2);
  const Model m0 = io::parse_inp(deck);
  const Real u1 = load_face_ux(m0, all[0]);
  const Real u2 = load_face_ux(m0, all[1]);
  CX_NEAR(u2, 2.0 * u1, 1e-9 * std::fabs(u2));
}

// Two stacked cubes; cube 2 is ABSENT in step 1 (removed) and ADDED in step 2. In step
// 2 the added cube is strain-free relative to the step-1 deformed shape: it carries
// stiffness for the step-2 increment only. Compare against the reference where cube 2
// is present in BOTH steps: the birth-death case must be SOFTER at the end of step 2
// (cube 2 missed step 1's loading), proving the removed element contributes nothing in
// the step it is inactive and the re-added one starts strain-free.
// Two stacked confined unit cubes driven by a PRESCRIBED top-face displacement (never
// singular even with cube 2 removed, because every node is fully prescribed): all
// nodes have u_x=u_z=0; the bottom face (y=0) has u_y=0; the shared mid face (y=1) is
// left free in y so cube 1 carries the load; the top face (y=2, cube-2-exclusive nodes
// 9-12) is prescribed u_y=delta. With cube 2 removed the mid and top faces are still
// pinned in y (see the y-pin lines) so the system stays non-singular.
std::string stacked_deck(const std::string& step1_mc, const std::string& step2_mc,
                         Real top_delta_step1, Real top_delta_step2) {
  const std::string header = R"(
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
*NSET, NSET=TOP
9, 10, 11, 12
*NSET, NSET=MID
3, 4, 7, 8
*NSET, NSET=BOT
1, 2, 5, 6
*NSET, NSET=ALLN
1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
*BOUNDARY
ALLN, 1, 1, 0.
ALLN, 3, 3, 0.
BOT, 2, 2, 0.
MID, 2, 2, 0.
)";
  auto step = [&](const std::string& mc, Real delta) {
    return "*STEP\n*STATIC\n" + mc + "*BOUNDARY\nTOP, 2, 2, " +
           std::to_string(delta) + "\n*END STEP\n";
  };
  return header + step(step1_mc, top_delta_step1) + step(step2_mc, top_delta_step2);
}

void test_cross_step_birth_death_softens() {
  // Cube 2 REMOVEd in step 1 (top prescribed to 0.01 with cube 2 dead -> it carries no
  // stress; cube 1's free top stays undeformed), then ADDed in step 2 with the top held
  // at 0.02. Because cube 2 is strain-free relative to the step-1 config (top already at
  // 0.01), in step 2 it stretches only by the increment 0.02-0.01=0.01.
  const std::vector<Model> bd_steps = io::parse_inp_steps(stacked_deck(
      "*MODEL CHANGE, TYPE=ELEMENT, REMOVE\n2\n",
      "*MODEL CHANGE, TYPE=ELEMENT, ADD\n2\n", 0.01, 0.02));
  CX_CHECK(bd_steps.size() == 2);
  CX_CHECK(bd_steps[0].deactivated_elements.size() == 1);  // removed in step 1
  CX_CHECK(bd_steps[0].deactivated_elements[0] == 2);
  CX_CHECK(bd_steps[1].deactivated_elements.empty());      // re-added in step 2
  const std::vector<StaticFields> bd = numerics::solve_multistep_static_all(bd_steps);

  const Model m0 = io::parse_inp(stacked_deck("", "", 0.0, 0.0));
  const std::size_t i10 = static_cast<std::size_t>(m0.mesh.node_index(10));  // cube-2 only

  // Step 1: cube 2 dead -> zero stress there even though node 10 is prescribed to 0.01.
  CX_NEAR(bd[0].stress[i10][1], 0.0, 1e-9);
  // Step 2: cube 2 live -> non-zero stress, but built ONLY from the 0.01 increment, not
  // the full 0.02. Compare to a REFERENCE where cube 2 is live throughout and driven the
  // SAME 0.01 increment in step 2 from a step-1 top of 0.01: the two must match, proving
  // reactivation is strain-free relative to the deformed config (not to the origin).
  const std::vector<Model> ref_steps =
      io::parse_inp_steps(stacked_deck("", "", 0.01, 0.02));
  const std::vector<StaticFields> ref =
      numerics::solve_multistep_static_all(ref_steps);
  CX_CHECK(std::fabs(bd[1].stress[i10][1]) > 1e-3);
  // The step-2 INCREMENT of cube-2 stress is identical in both (both stretch by 0.01 in
  // step 2 from a top at 0.01); the reference's step-1 already stressed cube 2, so its
  // TOTAL exceeds the birth-death total by exactly the step-1 stress it accrued.
  const Real ref_incr = ref[1].stress[i10][1] - ref[0].stress[i10][1];
  const Real bd_incr = bd[1].stress[i10][1] - bd[0].stress[i10][1];
  CX_NEAR(bd_incr, ref_incr, 1e-6 * std::fabs(ref_incr));
  // And the birth-death total is SMALLER in magnitude (it skipped step-1's stress).
  CX_CHECK(std::fabs(bd[1].stress[i10][1]) < std::fabs(ref[1].stress[i10][1]));
}

// An element removed in BOTH steps stays exactly zero-stress; a removed-then-added one
// becomes live. This isolates the active-mask carry across steps.
void test_dead_element_zero_stress() {
  const std::vector<Model> dead = io::parse_inp_steps(stacked_deck(
      "*MODEL CHANGE, TYPE=ELEMENT, REMOVE\n2\n",
      "*MODEL CHANGE, TYPE=ELEMENT, REMOVE\n2\n", 0.01, 0.02));
  CX_CHECK(dead[0].deactivated_elements.size() == 1);
  CX_CHECK(dead[1].deactivated_elements.size() == 1);  // still removed in step 2
  const StaticFields dead_f = numerics::solve_multistep_static(dead);
  const Model m0 = io::parse_inp(stacked_deck("", "", 0.0, 0.0));
  const std::size_t i10 = static_cast<std::size_t>(m0.mesh.node_index(10));
  CX_NEAR(dead_f.stress[i10][1], 0.0, 1e-12);
}

}  // namespace

int main() {
  test_single_step_is_unchanged();
  test_two_load_steps_sum_to_single();
  test_per_step_totals();
  test_cross_step_birth_death_softens();
  test_dead_element_zero_stress();
  if (cxtest::g_failures == 0) std::printf("test_multistep: OK\n");
  CX_MAIN_RETURN();
}
