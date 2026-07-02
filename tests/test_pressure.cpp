// *DLOAD pressure: consistent-load total (= -p*Area*n_hat) and solve-level global
// equilibrium (sum of reactions = -total applied load).
#include <array>
#include <cmath>

#include "calculixpp/fem/loads.hpp"
#include "calculixpp/numerics/linear_static.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

Model unit_tet(ElementType type) {
  Model m;
  m.mesh.add_node(1, {0, 0, 0});
  m.mesh.add_node(2, {1, 0, 0});
  m.mesh.add_node(3, {0, 1, 0});
  m.mesh.add_node(4, {0, 0, 1});
  std::vector<Index> conn = {1, 2, 3, 4};
  if (type == ElementType::C3D10) {
    // midside nodes 5..10 per nonei10 edges (1-2,2-3,3-1,1-4,2-4,3-4)
    m.mesh.add_node(5, {0.5, 0, 0});
    m.mesh.add_node(6, {0.5, 0.5, 0});
    m.mesh.add_node(7, {0, 0.5, 0});
    m.mesh.add_node(8, {0, 0, 0.5});
    m.mesh.add_node(9, {0.5, 0, 0.5});
    m.mesh.add_node(10, {0, 0.5, 0.5});
    conn = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  }
  m.mesh.add_element(100, type, conn);
  m.mesh.add_elset("EALL", {100});
  m.materials["EL"] = Material{"EL", ElasticIso{210000.0, 0.3}, std::nullopt};
  m.sections.push_back(SolidSection{"EALL", "EL"});
  return m;
}

// Face 1 (z=0 plane, area 1/2, outward normal -z): total force = (0,0,+p/2).
void test_consistent_total() {
  for (ElementType t : {ElementType::C3D4, ElementType::C3D10}) {
    Model m = unit_tet(t);
    m.dloads.push_back(Dload{100, 1, 100.0});
    const std::vector<Real> f = fem::external_load_vector(m);
    std::array<Real, 3> sum{0, 0, 0};
    for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i)
      for (int c = 0; c < 3; ++c) sum[static_cast<std::size_t>(c)] += f[i * 3 + static_cast<std::size_t>(c)];
    CX_NEAR(sum[0], 0.0, 1e-9);
    CX_NEAR(sum[1], 0.0, 1e-9);
    CX_NEAR(sum[2], 50.0, 1e-9);  // -p*Area*(-z) = +100*0.5
  }
}

// Fix nodes 1,2,3; pressure p=100 on face 3 (nodes 2,3,4, outward (1,1,1)/√3,
// area √3/2): total applied = (-50,-50,-50), so sum of reactions = (50,50,50).
void test_solve_equilibrium() {
  Model m = unit_tet(ElementType::C3D4);
  for (Index nd : {1, 2, 3})
    for (int c = 1; c <= 3; ++c) m.spcs.push_back(Spc{nd, c, 0.0});
  m.dloads.push_back(Dload{100, 3, 100.0});

  const StaticFields r = numerics::solve_linear_static(m);
  std::array<Real, 3> rf{0, 0, 0};
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i)
    for (int c = 0; c < 3; ++c) rf[static_cast<std::size_t>(c)] += r.reaction[i][static_cast<std::size_t>(c)];
  for (int c = 0; c < 3; ++c) CX_NEAR(rf[static_cast<std::size_t>(c)], 50.0, 1e-6);
  // Free node 4 should move under the pressure.
  const Real u4 = std::sqrt(r.displacement[3][0] * r.displacement[3][0] +
                            r.displacement[3][1] * r.displacement[3][1] +
                            r.displacement[3][2] * r.displacement[3][2]);
  CX_CHECK(u4 > 0.0);
}

}  // namespace

int main() {
  test_consistent_total();
  test_solve_equilibrium();
  if (cxtest::g_failures == 0) std::printf("test_pressure: OK\n");
  CX_MAIN_RETURN();
}
