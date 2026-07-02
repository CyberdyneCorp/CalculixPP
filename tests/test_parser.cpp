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

}  // namespace

int main() {
  test_embedded();
  test_surface();
  test_surface_errors();
  test_beam10p();
  if (cxtest::g_failures == 0) std::printf("test_parser: OK\n");
  CX_MAIN_RETURN();
}
