// *AMPLITUDE evaluation + amplitude-scaled loads, *DLOAD body loads (GRAV /
// CENTRIF), and *DSLOAD, exercising fem::external_load_vector and the parser.
#include <array>
#include <cmath>

#include "calculixpp/core/amplitude.hpp"
#include "calculixpp/fem/element.hpp"
#include "calculixpp/fem/loads.hpp"
#include "calculixpp/io/inp_parser.hpp"
#include "calculixpp/numerics/nonlinear_static.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

Model unit_tet() {
  Model m;
  m.mesh.add_node(1, {0, 0, 0});
  m.mesh.add_node(2, {1, 0, 0});
  m.mesh.add_node(3, {0, 1, 0});
  m.mesh.add_node(4, {0, 0, 1});
  m.mesh.add_element(100, ElementType::C3D4, {1, 2, 3, 4});
  m.mesh.add_elset("EALL", {100});
  m.materials["EL"] = Material{"EL", ElasticIso{210000.0, 0.3}, 7800.0};
  m.sections.push_back(SolidSection{"EALL", "EL"});
  return m;
}

// Tabular amplitude interpolates linearly between tabulated points.
void test_amplitude_tabular() {
  Amplitude a;
  a.definition = Amplitude::Definition::Tabular;
  a.points = {{0.0, 0.0}, {1.0, 2.0}, {2.0, 2.0}};
  CX_NEAR(a.value_at(0.0), 0.0, 1e-12);
  CX_NEAR(a.value_at(0.5), 1.0, 1e-12);   // halfway on the 0->2 ramp
  CX_NEAR(a.value_at(1.0), 2.0, 1e-12);
  CX_NEAR(a.value_at(1.5), 2.0, 1e-12);   // flat segment
  CX_NEAR(a.value_at(-1.0), 0.0, 1e-12);  // clamp below range
  CX_NEAR(a.value_at(5.0), 2.0, 1e-12);   // clamp above range
}

// STEP amplitude holds the left point's value.
void test_amplitude_step() {
  Amplitude a;
  a.definition = Amplitude::Definition::Step;
  a.points = {{0.0, 1.0}, {1.0, 3.0}};
  CX_NEAR(a.value_at(0.5), 1.0, 1e-12);  // still the first value
  CX_NEAR(a.value_at(1.0), 3.0, 1e-12);
}

// Periodic amplitude wraps time by its period.
void test_amplitude_periodic() {
  Amplitude a;
  a.definition = Amplitude::Definition::Tabular;
  a.points = {{0.0, 0.0}, {1.0, 1.0}};
  a.period = 1.0;
  CX_NEAR(a.value_at(0.25), 0.25, 1e-12);
  CX_NEAR(a.value_at(1.25), 0.25, 1e-12);  // wrapped into [0,1)
  CX_NEAR(a.value_at(2.75), 0.75, 1e-12);
}

// A *CLOAD with an amplitude is scaled by amplitude(step_fraction), not by the
// default linear ramp.
void test_amplitude_scaled_load() {
  Model m = unit_tet();
  Amplitude a;
  a.points = {{0.0, 0.0}, {1.0, 1.0}, {2.0, 1.0}};
  a.name = "RAMP";
  m.amplitudes["RAMP"] = a;
  m.increment.total = 2.0;                        // step spans t in [0,2]
  m.cloads.push_back(Cload{4, 3, 1000.0, "RAMP"});  // Fz=1000 on node 4

  // At lambda=0.25 -> t=0.5 -> amplitude 0.5 -> applied 500.
  const std::vector<Real> f = fem::external_load_vector(m, 0.25);
  CX_NEAR(f[3 * 3 + 2], 500.0, 1e-9);
  // Full step (lambda=1 -> t=2 -> amplitude 1.0) -> full 1000.
  const std::vector<Real> ff = fem::external_load_vector(m);
  CX_NEAR(ff[3 * 3 + 2], 1000.0, 1e-9);

  // A load WITHOUT an amplitude keeps the linear ramp: half at lambda=0.5.
  Model m2 = unit_tet();
  m2.cloads.push_back(Cload{4, 3, 1000.0, ""});
  const std::vector<Real> g = fem::external_load_vector(m2, 0.5);
  CX_NEAR(g[3 * 3 + 2], 500.0, 1e-9);
}

// Gravity on a fixed block: total reaction magnitude = rho * Volume * g.
void test_gravity_equilibrium() {
  Model m = unit_tet();
  const Real g = 9.81, rho = 7800.0;
  BodyLoad bl;
  bl.kind = BodyLoad::Kind::Gravity;
  bl.magnitude = g;
  bl.dir = {0, 0, -1};
  m.body_loads.push_back(bl);

  const Real vol = fem::element_volume(ElementType::C3D4,
                                       std::array<Vec3, 4>{
                                           Vec3{0, 0, 0}, Vec3{1, 0, 0},
                                           Vec3{0, 1, 0}, Vec3{0, 0, 1}});
  const std::vector<Real> f = fem::external_load_vector(m);
  Real sumz = 0.0;
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i) sumz += f[i * 3 + 2];
  CX_NEAR(sumz, -rho * vol * g, 1e-6);  // total gravity force, downward

  // Fix all nodes and solve: sum of reactions balances the applied gravity.
  for (Index nd : {1, 2, 3, 4})
    for (int c = 1; c <= 3; ++c) m.spcs.push_back(Spc{nd, c, 0.0});
  const StaticFields r = numerics::solve_nonlinear_static(m);
  Real rz = 0.0;
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i) rz += r.reaction[i][2];
  CX_NEAR(rz, rho * vol * g, 1e-4);  // reaction = -applied
}

// Centrifugal load: finite and radially outward from the axis.
void test_centrifugal_radial() {
  Model m = unit_tet();
  BodyLoad bl;
  bl.kind = BodyLoad::Kind::Centrifugal;
  bl.magnitude = 100.0;      // omega^2
  bl.point = {0, 0, 0};      // axis through origin
  bl.dir = {0, 0, 1};        // spin about z -> radial in x,y
  m.body_loads.push_back(bl);

  const std::vector<Real> f = fem::external_load_vector(m);
  // The resultant in-plane force must point away from the z-axis (outward): the
  // centroid is in the +x,+y octant, so sum Fx, Fy > 0, and Fz ~ 0.
  Real fx = 0, fy = 0, fz = 0;
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i) {
    fx += f[i * 3 + 0];
    fy += f[i * 3 + 1];
    fz += f[i * 3 + 2];
  }
  CX_CHECK(std::isfinite(fx) && std::isfinite(fy));
  CX_CHECK(fx > 0.0);
  CX_CHECK(fy > 0.0);
  CX_NEAR(fz, 0.0, 1e-9);  // no axial component
}

// *DSLOAD reuses the pressure-face machinery: same nodal forces as *DLOAD P<face>.
void test_dsload_matches_dload() {
  const char* deck_d =
      "*NODE\n1,0,0,0\n2,1,0,0\n3,0,1,0\n4,0,0,1\n"
      "*ELEMENT,TYPE=C3D4,ELSET=EALL\n100,1,2,3,4\n"
      "*MATERIAL,NAME=EL\n*ELASTIC\n210000,0.3\n"
      "*SOLID SECTION,ELSET=EALL,MATERIAL=EL\n"
      "*STEP\n*STATIC\n*DLOAD\n100,P1,100.0\n*END STEP\n";
  const char* deck_s =
      "*NODE\n1,0,0,0\n2,1,0,0\n3,0,1,0\n4,0,0,1\n"
      "*ELEMENT,TYPE=C3D4,ELSET=EALL\n100,1,2,3,4\n"
      "*MATERIAL,NAME=EL\n*ELASTIC\n210000,0.3\n"
      "*SOLID SECTION,ELSET=EALL,MATERIAL=EL\n"
      "*STEP\n*STATIC\n*DSLOAD\n100,P1,100.0\n*END STEP\n";
  const Model md = io::parse_inp(deck_d);
  const Model ms = io::parse_inp(deck_s);
  const std::vector<Real> fd = fem::external_load_vector(md);
  const std::vector<Real> fs = fem::external_load_vector(ms);
  CX_CHECK(fd.size() == fs.size());
  for (std::size_t i = 0; i < fd.size(); ++i) CX_NEAR(fd[i], fs[i], 1e-12);
}

// Parser wires *AMPLITUDE, AMPLITUDE= on *CLOAD, and *DLOAD GRAV end to end.
void test_parse_amplitude_and_gravity() {
  const char* deck =
      "*NODE\n1,0,0,0\n2,1,0,0\n3,0,1,0\n4,0,0,1\n"
      "*ELEMENT,TYPE=C3D4,ELSET=EALL\n100,1,2,3,4\n"
      "*MATERIAL,NAME=EL\n*ELASTIC\n210000,0.3\n*DENSITY\n7800\n"
      "*SOLID SECTION,ELSET=EALL,MATERIAL=EL\n"
      "*AMPLITUDE,NAME=A1,DEFINITION=TABULAR\n0.0,0.0,1.0,1.0\n"
      "*STEP\n*STATIC\n"
      "*CLOAD,AMPLITUDE=A1\n4,3,1000.0\n"
      "*DLOAD\nEALL,GRAV,9.81,0,0,-1\n*END STEP\n";
  const Model m = io::parse_inp(deck);
  CX_CHECK(m.amplitudes.count("A1") == 1);
  CX_CHECK(m.amplitudes.at("A1").points.size() == 2);
  CX_CHECK(m.cloads.size() == 1);
  CX_CHECK(m.cloads[0].amplitude == "A1");
  CX_CHECK(m.body_loads.size() == 1);
  CX_CHECK(m.body_loads[0].kind == BodyLoad::Kind::Gravity);
  CX_NEAR(m.body_loads[0].magnitude, 9.81, 1e-12);

  // Half step: cload amplitude(0.5)=0.5 -> 500; gravity keeps the linear ramp 0.5.
  // Total Fz = 500 (cload) - 0.5 * rho*V*g (gravity), checked as the global sum.
  const std::vector<Real> f = fem::external_load_vector(m, 0.5);
  const Real vol = 1.0 / 6.0;  // unit tet volume
  Real fz = 0.0;
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i) fz += f[i * 3 + 2];
  CX_NEAR(fz, 500.0 - 0.5 * 7800.0 * vol * 9.81, 1e-4);
}

}  // namespace

int main() {
  test_amplitude_tabular();
  test_amplitude_step();
  test_amplitude_periodic();
  test_amplitude_scaled_load();
  test_gravity_equilibrium();
  test_centrifugal_radial();
  test_dsload_matches_dload();
  test_parse_amplitude_and_gravity();
  if (cxtest::g_failures == 0) std::printf("test_loads: OK\n");
  CX_MAIN_RETURN();
}
