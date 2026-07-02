// Load-accumulation (OP=MOD / OP=NEW) and step property-change cards
// (*CHANGE MATERIAL / *CHANGE PLASTIC / *CHANGE SOLID SECTION), parsed and
// applied within the single-step model. (tasks.md 2.3 / 2.4.)
#include <string>

#include "calculixpp/io/inp_parser.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

// A minimal deck header shared by the load tests: one tet, a material and a
// solid section, with named sets so loads can be expressed compactly.
const char* kHeader = R"(
*NODE, NSET=NALL
1, 0., 0., 0.
2, 1., 0., 0.
3, 0., 1., 0.
4, 0., 0., 1.
*ELEMENT, TYPE=C3D4, ELSET=EALL
100, 1, 2, 3, 4
*MATERIAL, NAME=EL
*ELASTIC
210000., 0.3
*DENSITY
7800.
*SOLID SECTION, ELSET=EALL, MATERIAL=EL
*STEP
*STATIC
)";

Model parse(const std::string& body) {
  return io::parse_inp(std::string(kHeader) + body + "\n*END STEP\n");
}

// OP defaults to MOD: successive *CLOAD cards accumulate (both survive).
void test_cload_mod_accumulates() {
  const Model m = parse(R"(
*CLOAD
4, 3, 10.
*CLOAD
4, 2, 5.
)");
  CX_CHECK(m.cloads.size() == 2);
}

// OP=NEW on a later *CLOAD resets the concentrated loads: the earlier one is
// dropped, only the OP=NEW card's load remains.
void test_cload_new_resets() {
  const Model m = parse(R"(
*CLOAD
4, 3, 10.
*CLOAD, OP=NEW
4, 2, 5.
)");
  CX_CHECK(m.cloads.size() == 1);
  CX_CHECK(m.cloads[0].comp == 2);
  CX_NEAR(m.cloads[0].value, 5.0, 1e-12);
}

// OP=NEW takes effect on the FIRST card of the type in the step only: a second
// OP=NEW behaves like MOD and does not wipe what the first one established.
void test_cload_new_first_card_only() {
  const Model m = parse(R"(
*CLOAD, OP=NEW
4, 3, 10.
*CLOAD, OP=NEW
4, 2, 5.
)");
  CX_CHECK(m.cloads.size() == 2);
}

// A multi-line *CLOAD card with OP=NEW keeps all of its own lines (the reset
// fires once per card at the header, not per data line).
void test_cload_new_multiline_card() {
  const Model m = parse(R"(
*CLOAD
2, 1, 1.
*CLOAD, OP=NEW
4, 3, 10.
4, 2, 5.
)");
  CX_CHECK(m.cloads.size() == 2);  // the pre-NEW load dropped, both NEW lines kept
}

// OP=NEW on *BOUNDARY resets prescribed DOFs.
void test_boundary_new_resets() {
  const Model m = parse(R"(
*BOUNDARY
1, 1, 3
*BOUNDARY, OP=NEW
2, 1, 1
)");
  CX_CHECK(m.spcs.size() == 1);  // 3 DOFs of node 1 dropped, one DOF of node 2 kept
  CX_CHECK(m.spcs[0].node_id == 2);
}

// OP=NEW on *DLOAD clears both pressure and body loads of the distributed set.
void test_dload_new_resets() {
  const Model m = parse(R"(
*DLOAD
100, P1, 2.
*DLOAD
EALL, GRAV, 9.81, 0., 0., -1.
*DLOAD, OP=NEW
100, P2, 3.
)");
  CX_CHECK(m.dloads.size() == 1);       // P1 pressure dropped, P2 kept
  CX_CHECK(m.dloads[0].face == 2);
  CX_CHECK(m.body_loads.empty());       // gravity dropped by OP=NEW
}

// An unknown OP value is rejected.
void test_op_invalid_throws() {
  bool threw = false;
  try {
    parse("*CLOAD, OP=FOO\n4, 3, 10.\n");
  } catch (const io::ParseError&) {
    threw = true;
  }
  CX_CHECK(threw);
}

// *CHANGE SOLID SECTION re-binds a material to an element set within the step;
// the later binding wins (element_elastic is last-writer per element).
void test_change_solid_section_rebinds() {
  // SOFT is declared before the step; the *CHANGE SOLID SECTION inside the step
  // re-binds EALL to it. Material declaration order is independent of section
  // resolution, which happens in element_elastic().
  const Model m = io::parse_inp(
      "*MATERIAL, NAME=SOFT\n*ELASTIC\n1000., 0.3\n" + std::string(kHeader) +
      "*CHANGE SOLID SECTION, ELSET=EALL, MATERIAL=SOFT\n*END STEP\n");
  const auto el = m.element_elastic();
  CX_CHECK(el.size() == 1);
  CX_NEAR(el[0].E, 1000.0, 1e-9);  // SOFT overrides the header's EL (210000)
}

// *CHANGE SOLID SECTION requires MATERIAL= and ELSET=.
void test_change_solid_section_requires_params() {
  bool threw = false;
  try {
    parse("*CHANGE SOLID SECTION, ELSET=EALL\n");
  } catch (const io::ParseError&) {
    threw = true;
  }
  CX_CHECK(threw);
}

// *CHANGE MATERIAL re-opens an existing material; unknown names are rejected.
// The following *CHANGE PLASTIC data is accepted and ignored (no plasticity
// yet) — the deck parses without error.
void test_change_material_plastic_parses() {
  const Model m = parse(R"(
*CHANGE MATERIAL, NAME=EL
*CHANGE PLASTIC
0., 0.
0.01, 100.
)");
  // Elastic data is unchanged (plastic change is a no-op for now).
  const auto el = m.element_elastic();
  CX_NEAR(el[0].E, 210000.0, 1e-6);
}

void test_change_material_unknown_throws() {
  bool threw = false;
  try {
    parse("*CHANGE MATERIAL, NAME=NOPE\n*CHANGE PLASTIC\n0.,0.\n");
  } catch (const io::ParseError&) {
    threw = true;
  }
  CX_CHECK(threw);
}

}  // namespace

int main() {
  test_cload_mod_accumulates();
  test_cload_new_resets();
  test_cload_new_first_card_only();
  test_cload_new_multiline_card();
  test_boundary_new_resets();
  test_dload_new_resets();
  test_op_invalid_throws();
  test_change_solid_section_rebinds();
  test_change_solid_section_requires_params();
  test_change_material_plastic_parses();
  test_change_material_unknown_throws();
  if (cxtest::g_failures == 0) std::printf("test_step_changes: OK\n");
  CX_MAIN_RETURN();
}
