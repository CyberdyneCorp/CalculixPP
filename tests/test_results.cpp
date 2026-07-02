// Stress recovery (uniaxial-stress patch test) and .frd/.dat writer smoke test.
#include <fstream>
#include <sstream>
#include <string>

#include "calculixpp/fem/stress.hpp"
#include "calculixpp/io/results_writer.hpp"
#include "calculixpp/numerics/linear_static.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

Model tet_model(bool loaded) {
  Model m;
  m.mesh.add_node(1, {0, 0, 0});
  m.mesh.add_node(2, {2, 0, 0});
  m.mesh.add_node(3, {0.3, 1.5, 0});
  m.mesh.add_node(4, {0.1, 0.2, 3});
  m.mesh.add_element(100, ElementType::C3D4, {1, 2, 3, 4});
  m.mesh.add_elset("EALL", {100});
  m.materials["EL"] = Material{"EL", ElasticIso{210000.0, 0.3}, std::nullopt};
  m.sections.push_back(SolidSection{"EALL", "EL"});
  if (loaded) {
    for (Index nd : {1, 2, 3})
      for (int c = 1; c <= 3; ++c) m.spcs.push_back(Spc{nd, c, 0.0});
    m.cloads.push_back(Cload{4, 3, -1000.0});
  }
  return m;
}

// Impose a uniaxial-stress displacement field; recovered stress must be uniform.
void test_uniaxial_patch() {
  const Model m = tet_model(false);
  const Real E = 210000.0, nu = 0.3, S = 100.0;
  const Real exx = S / E, elat = -nu * S / E;

  StaticFields f;
  f.displacement.resize(m.mesh.num_nodes());
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i) {
    const Vec3 x = m.mesh.nodes()[i].x;
    f.displacement[i] = {exx * x[0], elat * x[1], elat * x[2]};
  }
  fem::recover_fields(m, f);

  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i) {
    CX_NEAR(f.stress[i][0], S, 1e-6);      // SXX = S
    for (int c = 1; c < 6; ++c) CX_NEAR(f.stress[i][static_cast<std::size_t>(c)], 0.0, 1e-6);
  }
}

void test_writers() {
  const Model m = tet_model(true);
  const StaticFields f = numerics::solve_linear_static(m);

  io::write_frd("cxpp_out.frd", m, f);
  io::write_dat("cxpp_out.dat", m, f);

  std::ifstream frd("cxpp_out.frd"), dat("cxpp_out.dat");
  CX_CHECK(frd.good());
  CX_CHECK(dat.good());
  std::stringstream fs, ds;
  fs << frd.rdbuf();
  ds << dat.rdbuf();
  CX_CHECK(fs.str().find(" -3") != std::string::npos);       // frd block terminator
  CX_CHECK(fs.str().find("DISP") != std::string::npos);
  CX_CHECK(fs.str().find("STRESS") != std::string::npos);
  CX_CHECK(ds.str().find("displacements") != std::string::npos);
}

}  // namespace

int main() {
  test_uniaxial_patch();
  test_writers();
  if (cxtest::g_failures == 0) std::printf("test_results: OK\n");
  CX_MAIN_RETURN();
}
