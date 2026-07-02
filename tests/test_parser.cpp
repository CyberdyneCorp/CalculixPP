// Parser tests: an embedded single-tet deck (structure + solve) and, when
// available, the real beam10p.inp cantilever parsed and solved end to end.
#include <fstream>
#include <string>

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

  const numerics::LinearStaticResult r = numerics::solve_linear_static(m);
  CX_CHECK(r.displacement.size() == 4);
  CX_CHECK(r.displacement[3][2] < 0.0);  // apex pushed in -z
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

  const numerics::LinearStaticResult r = numerics::solve_linear_static(m);
  // Cantilever loaded in +y at the free (z=8) end -> tip deflects +y.
  const Index tip = m.mesh.node_index(11);
  CX_CHECK(tip >= 0);
  CX_CHECK(r.displacement[static_cast<std::size_t>(tip)][1] > 0.0);
}

}  // namespace

int main() {
  test_embedded();
  test_beam10p();
  if (cxtest::g_failures == 0) std::printf("test_parser: OK\n");
  CX_MAIN_RETURN();
}
