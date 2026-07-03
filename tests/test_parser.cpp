// Parser tests: an embedded single-tet deck (structure + solve) and, when
// available, the real beam10p.inp cantilever parsed and solved end to end.
#include <fstream>
#include <string>
#include <utility>

#include "calculixpp/io/inp_parser.hpp"
#include "calculixpp/numerics/linear_static.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

const char* kSingleTet = R"(
** single C3D4 tet, fixed base, point load on the apex
*NODE, NSET=NALL
1, 0., 0., 0.
2, 1., 0., 0.
3, 0., 1., 0.
4, 0., 0., 1.
*ELEMENT, TYPE=C3D4, ELSET=EALL
100, 1, 2, 3, 4
*NSET, NSET=FIX
1, 2, 3
*BOUNDARY
FIX, 1, 3
*MATERIAL, NAME=EL
*ELASTIC
210000., 0.3
*SOLID SECTION, ELSET=EALL, MATERIAL=EL
*STEP
*STATIC
*CLOAD
4, 3, -1000.
*END STEP
)";

void test_embedded() {
  const Model m = io::parse_inp(kSingleTet);
  CX_CHECK(m.mesh.num_nodes() == 4);
  CX_CHECK(m.mesh.num_elements() == 1);
  CX_CHECK(m.mesh.elements()[0].type == ElementType::C3D4);
  CX_CHECK(m.materials.count("EL") == 1);
  CX_CHECK(m.materials.at("EL").elastic.has_value());
  CX_NEAR(m.materials.at("EL").elastic->E, 210000.0, 1e-9);
  CX_NEAR(m.materials.at("EL").elastic->nu, 0.3, 1e-12);
  CX_CHECK(m.sections.size() == 1);
  CX_CHECK(m.spcs.size() == 9);   // 3 nodes x 3 dofs
  CX_CHECK(m.cloads.size() == 1);
  CX_CHECK(m.mesh.nset("FIX") != nullptr);

  const StaticFields r = numerics::solve_linear_static(m);
  CX_CHECK(r.displacement.size() == 4);
  CX_CHECK(r.displacement[3][2] < 0.0);  // apex pushed in -z
}

const char* kSurfaceDeck = R"(
*NODE, NSET=NALL
1, 0., 0., 0.
2, 1., 0., 0.
3, 0., 1., 0.
4, 0., 0., 1.
5, 1., 1., 1.
*ELEMENT, TYPE=C3D4, ELSET=EALL
100, 1, 2, 3, 4
101, 2, 3, 4, 5
*ELSET, ELSET=LOADED
100, 101
*NSET, NSET=TOP
4, 5
*SURFACE, NAME=SPRESS, TYPE=ELEMENT
LOADED, S2
100, S3
*SURFACE, NAME=SDEFAULT
101, S1
*SURFACE, NAME=SNODES, TYPE=NODE
TOP
1
)";

void test_surface() {
  const Model m = io::parse_inp(kSurfaceDeck);

  // Element surface expanding an elset ref plus an explicit element.
  const Surface* sp = m.mesh.surface("SPRESS");
  CX_CHECK(sp != nullptr);
  CX_CHECK(sp->type == Surface::Type::Element);
  CX_CHECK(sp->faces.size() == 3);  // 100/S2, 101/S2 (from elset), 100/S3
  CX_CHECK((sp->faces[0] == std::pair<Index, int>{100, 2}));
  CX_CHECK((sp->faces[1] == std::pair<Index, int>{101, 2}));
  CX_CHECK((sp->faces[2] == std::pair<Index, int>{100, 3}));

  // TYPE defaults to ELEMENT.
  const Surface* sd = m.mesh.surface("SDEFAULT");
  CX_CHECK(sd != nullptr);
  CX_CHECK(sd->type == Surface::Type::Element);
  CX_CHECK(sd->faces.size() == 1);
  CX_CHECK((sd->faces[0] == std::pair<Index, int>{101, 1}));

  // Node surface expanding an nset ref plus an explicit node id.
  const Surface* sn = m.mesh.surface("SNODES");
  CX_CHECK(sn != nullptr);
  CX_CHECK(sn->type == Surface::Type::Node);
  CX_CHECK(sn->nodes.size() == 3);  // 4, 5 (from nset TOP), 1
  CX_CHECK(sn->nodes[0] == 4);
  CX_CHECK(sn->nodes[1] == 5);
  CX_CHECK(sn->nodes[2] == 1);

  CX_CHECK(m.mesh.surface("MISSING") == nullptr);
  CX_CHECK(m.mesh.surfaces().size() == 3);
}

void test_surface_errors() {
  bool threw = false;
  try {
    io::parse_inp("*SURFACE, NAME=BAD\n100, X7\n");
  } catch (const io::ParseError&) {
    threw = true;
  }
  CX_CHECK(threw);  // face label must be S<n>

  threw = false;
  try {
    io::parse_inp("*SURFACE, NAME=BAD, TYPE=FOO\n");
  } catch (const io::ParseError&) {
    threw = true;
  }
  CX_CHECK(threw);  // unsupported TYPE
}

void test_beam10p() {
  const std::string path = "/home/leonardo/work/CalculiX/test/beam10p.inp";
  std::ifstream probe(path);
  if (!probe) {
    std::printf("test_parser: beam10p.inp not present, skipping reference deck\n");
    return;
  }
  const Model m = io::parse_inp_file(path);
  CX_CHECK(m.mesh.num_nodes() == 90);
  CX_CHECK(m.mesh.num_elements() == 31);

  const StaticFields r = numerics::solve_linear_static(m);
  // Cantilever loaded in +y at the free (z=8) end -> tip deflects +y.
  const Index tip = m.mesh.node_index(11);
  CX_CHECK(tip >= 0);
  CX_CHECK(r.displacement[static_cast<std::size_t>(tip)][1] > 0.0);
}

// *PLASTIC / *CYCLIC HARDENING parsing: the hardening table lands on the material,
// HARDENING= selects the kind, and a plastic material flips has_plasticity().
void test_plastic_card() {
  const char* deck = R"(
*NODE
1, 0., 0., 0.
2, 1., 0., 0.
3, 0., 1., 0.
4, 0., 0., 1.
*ELEMENT, TYPE=C3D4, ELSET=EALL
100, 1, 2, 3, 4
*MATERIAL, NAME=STEEL
*ELASTIC
210000., 0.3
*PLASTIC, HARDENING=ISOTROPIC
800., 0.
960., 0.02
*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL
)";
  const Model m = io::parse_inp(deck);
  CX_CHECK(m.has_plasticity());
  const auto& mat = m.materials.at("STEEL");
  CX_CHECK(mat.plastic.has_value());
  CX_CHECK(mat.plastic->hardening == Plastic::Hardening::Isotropic);
  CX_CHECK(mat.plastic->yield.size() == 2);
  CX_NEAR(mat.plastic->yield[0], 800.0, 1e-12);
  CX_NEAR(mat.plastic->eqplastic[1], 0.02, 1e-12);

  // element_plastic() aligns with mesh.elements(); the single element is plastic.
  const auto ep = m.element_plastic();
  CX_CHECK(ep.size() == 1);
  CX_CHECK(ep[0].has_value());

  // An elastic-only deck stays non-plastic (linear path unaffected).
  const char* elastic = R"(
*NODE
1, 0., 0., 0.
*MATERIAL, NAME=EL
*ELASTIC
210000., 0.3
)";
  CX_CHECK(!io::parse_inp(elastic).has_plasticity());
}

// *USER MATERIAL / *DEPVAR / *RATEDEPENDENT and *HYPERELASTIC parsing (spec 4.3/4.6).
void test_user_and_hyperelastic_cards() {
  const char* user = R"(
*NODE
1, 0., 0., 0.
2, 1., 0., 0.
3, 0., 1., 0.
4, 0., 0., 1.
*ELEMENT, TYPE=C3D4, ELSET=EALL
100, 1, 2, 3, 4
*MATERIAL, NAME=UMAT
*USER MATERIAL, CONSTANTS=2, NAME=LINEAR_ELASTIC_USER
210000., 0.3
*DEPVAR
3
*RATEDEPENDENT
1.5
*SOLID SECTION, ELSET=EALL, MATERIAL=UMAT
)";
  const Model mu = io::parse_inp(user);
  CX_CHECK(mu.has_nonlinear_material());
  const auto& um = mu.materials.at("UMAT").user;
  CX_CHECK(um.has_value());
  CX_CHECK(um->name == "LINEAR_ELASTIC_USER");
  CX_CHECK(um->constants.size() == 2);
  CX_NEAR(um->constants[0], 210000.0, 1e-9);
  CX_CHECK(um->ndepvar == 3);
  CX_CHECK(um->rate_scale.has_value());
  CX_NEAR(*um->rate_scale, 1.5, 1e-12);
  // effective_elastic() derives (E,nu) from the user constants (no *ELASTIC given).
  const auto eff = mu.materials.at("UMAT").effective_elastic();
  CX_CHECK(eff.has_value());
  CX_NEAR(eff->E, 210000.0, 1e-6);

  const char* hyper = R"(
*NODE
1, 0., 0., 0.
2, 1., 0., 0.
3, 0., 1., 0.
4, 0., 0., 1.
*ELEMENT, TYPE=C3D4, ELSET=EALL
100, 1, 2, 3, 4
*MATERIAL, NAME=RUBBER
*HYPERELASTIC, N=1
0.5, 0.001
*SOLID SECTION, ELSET=EALL, MATERIAL=RUBBER
)";
  const Model mh = io::parse_inp(hyper);
  CX_CHECK(mh.has_nonlinear_material());
  const auto& hy = mh.materials.at("RUBBER").hyperelastic;
  CX_CHECK(hy.has_value());
  CX_NEAR(hy->c10, 0.5, 1e-12);
  CX_NEAR(hy->mu, 1.0, 1e-12);       // mu = 2*C10
  CX_NEAR(hy->kappa, 2000.0, 1e-9);  // kappa = 2/D1
}

// Parse the constraint cards and check they populate the Model and expand correctly.
void test_constraint_cards() {
  const char* deck = R"(
*NODE
1, 0., 0., 0.
2, 1., 0., 0.
3, 1., 1., 0.
4, 0., 1., 0.
5, 0., 0., 1.
6, 1., 0., 1.
7, 1., 1., 1.
8, 0., 1., 1.
*ELEMENT, TYPE=C3D8, ELSET=EALL
1, 1,2,3,4,5,6,7,8
*NSET, NSET=TOP
5,6,7,8
*EQUATION
2
6, 1, 1.0, 7, 1, -1.0
*MPC
BEAM, 5, 6
*RIGID BODY, NSET=TOP, REF NODE=5
*COUPLING, REF NODE=5, SURFACE=TOP
*KINEMATIC
1,3
*MATERIAL, NAME=STEEL
*ELASTIC
210000., 0.3
*SOLID SECTION, MATERIAL=STEEL, ELSET=EALL
)";
  const Model m = io::parse_inp(deck);
  CX_CHECK(m.equations.size() == 1);
  CX_CHECK(m.equations[0].terms.size() == 2);
  CX_CHECK(m.equations[0].terms[0].node_id == 6);
  CX_NEAR(m.equations[0].terms[0].coeff, 1.0, 1e-12);
  CX_CHECK(m.mpcs.size() == 1);
  CX_CHECK(m.mpcs[0].kind == Mpc::Kind::Beam);
  CX_CHECK(m.rigid_bodies.size() == 1);
  CX_CHECK(m.rigid_bodies[0].ref_node == 5);
  CX_CHECK(m.couplings.size() == 1);
  CX_CHECK(m.couplings[0].kind == Coupling::Kind::Kinematic);
  CX_CHECK(m.couplings[0].ref_node == 5);
  // Expansion turns the four constraint kinds into a flat equation list.
  const std::vector<Equation> eqs = m.expand_constraints();
  CX_CHECK(eqs.size() >= 1 + 3 + 3 + 9);  // eq + beam(3) + rigid(3 nodes*3) + coupling(3*3)
}

}  // namespace

int main() {
  test_embedded();
  test_surface();
  test_surface_errors();
  test_beam10p();
  test_plastic_card();
  test_user_and_hyperelastic_cards();
  test_constraint_cards();
  if (cxtest::g_failures == 0) std::printf("test_parser: OK\n");
  CX_MAIN_RETURN();
}
