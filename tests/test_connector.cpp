// Discrete/connector elements (*SPRING/*MASS/*DASHPOT). Validates the spring
// stiffness contribution against the analytical single-spring solution u = F/k for
// all three spring kinds, and that *MASS/*DASHPOT properties are parsed and stored
// (they are inert in a static solve, consumed by dynamics later).
#include <cmath>
#include <string>

#include "calculixpp/io/inp_parser.hpp"
#include "calculixpp/numerics/linear_static.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

const std::string kSpringA = R"(
*NODE,NSET=NALL
1,0.,0.,0.
2,1.,0.,0.
*ELEMENT,TYPE=SPRINGA,ELSET=EALL
1,1,2
*BOUNDARY
1,1,3
2,2,3
*SPRING,ELSET=EALL

10.
*STEP
*STATIC
*CLOAD
2,1,1.
*END STEP
)";

const std::string kSpring1 = R"(
*NODE,NSET=NALL
2,1.,0.,0.
*ELEMENT,TYPE=SPRING1,ELSET=EALL
1,2
*BOUNDARY
2,2,3
*SPRING,ELSET=EALL
1
10.
*STEP
*STATIC
*CLOAD
2,1,1.
*END STEP
)";

const std::string kSpring2 = R"(
*NODE,NSET=NALL
1,0.,0.,0.
2,1.,0.,0.
*ELEMENT,TYPE=SPRING2,ELSET=EALL
1,1,2
*BOUNDARY
1,1,3
2,2,3
*SPRING,ELSET=EALL
1,1
10.
*STEP
*STATIC
*CLOAD
2,1,2.
*END STEP
)";

// A model with a point mass and a dashpot alongside a spring — the static solve
// must ignore mass/dashpot but still store them for dynamics.
const std::string kMassDashpot = R"(
*NODE,NSET=NALL
1,0.,0.,0.
2,1.,0.,0.
*ELEMENT,TYPE=SPRINGA,ELSET=ESPRING
1,1,2
*ELEMENT,TYPE=MASS,ELSET=EMASS
2,2
*ELEMENT,TYPE=DASHPOTA,ELSET=EDASH
3,1,2
*SPRING,ELSET=ESPRING

10.
*MASS,ELSET=EMASS
0.9
*DASHPOT,ELSET=EDASH

1.e-5
*BOUNDARY
1,1,3
2,2,3
*STEP
*STATIC
*CLOAD
2,1,1.
*END STEP
)";

Real ux(const StaticFields& r, const Model& m, Index node) {
  return r.displacement[static_cast<std::size_t>(m.mesh.node_index(node))][0];
}

void test_spring_kinds() {
  const Model ma = io::parse_inp(kSpringA);
  CX_NEAR(ux(numerics::solve_linear_static(ma), ma, 2), 0.1, 1e-12);

  const Model m1 = io::parse_inp(kSpring1);
  CX_NEAR(ux(numerics::solve_linear_static(m1), m1, 2), 0.1, 1e-12);

  const Model m2 = io::parse_inp(kSpring2);
  CX_NEAR(ux(numerics::solve_linear_static(m2), m2, 2), 0.2, 1e-12);
}

void test_mass_dashpot_stored() {
  const Model m = io::parse_inp(kMassDashpot);
  // Spring still solves (u = 0.1); mass/dashpot are stored but inert statically.
  CX_NEAR(ux(numerics::solve_linear_static(m), m, 2), 0.1, 1e-12);
  CX_CHECK(m.springs.size() == 1);
  CX_CHECK(m.point_masses.size() == 1);
  CX_NEAR(m.point_masses[0].mass, 0.9, 1e-12);
  CX_CHECK(m.point_masses[0].node == 2);
  CX_CHECK(m.dashpots.size() == 1);
  CX_NEAR(m.dashpots[0].coefficient, 1.e-5, 1e-15);
}

}  // namespace

int main() {
  test_spring_kinds();
  test_mass_dashpot_stored();
  if (cxtest::g_failures == 0) std::printf("test_connector: OK\n");
  CX_MAIN_RETURN();
}
