// Steady-state heat conduction (spec: heat-transfer-analysis, task 3.1).
//   * conduction element matrix properties (row sums vanish: a constant temperature
//     field carries no flux);
//   * analytical 1-D bar: prescribed end temperatures -> linear profile + the exact
//     Fourier heat rate Q = k A dT / L recovered as the nodal flux reaction;
//   * thermal load vectors (*CFLUX nodal, *DFLUX consistent surface integral);
//   * end-to-end solve of a two-element bar reproduces the linear profile;
//   * the reference single-element decks (oneel20cf / oneel20df) are matched by the
//     Python regression suite to machine precision (see python/tests).
#include <array>
#include <cmath>
#include <numeric>
#include <vector>

#include "calculixpp/fem/element.hpp"
#include "calculixpp/fem/loads.hpp"
#include "calculixpp/io/inp_parser.hpp"
#include "calculixpp/numerics/heat_transfer.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

// A unit cube meshed as one C3D8 hex, corners in CalculiX order (bottom CCW then
// top CCW), with x the conduction axis. k is the isotropic conductivity.
Model unit_cube_hex(Real k) {
  Model m;
  // bottom face (z=0): 1..4, top face (z=1): 5..8. CalculiX C3D8 ordering.
  m.mesh.add_node(1, {0, 0, 0});
  m.mesh.add_node(2, {1, 0, 0});
  m.mesh.add_node(3, {1, 1, 0});
  m.mesh.add_node(4, {0, 1, 0});
  m.mesh.add_node(5, {0, 0, 1});
  m.mesh.add_node(6, {1, 0, 1});
  m.mesh.add_node(7, {1, 1, 1});
  m.mesh.add_node(8, {0, 1, 1});
  m.mesh.add_element(1, ElementType::C3D8, {1, 2, 3, 4, 5, 6, 7, 8});
  m.mesh.add_elset("EALL", {1});
  Material mat;
  mat.name = "EL";
  mat.thermal = Thermal{k, 0.0};
  m.materials["EL"] = mat;
  m.sections.push_back(SolidSection{"EALL", "EL"});
  m.procedure = Procedure::HeatTransferSteady;
  return m;
}

// The conduction matrix of any element has zero row sums: Σ_b Kt[a][b] = 0, because
// a uniform temperature (all nodal T equal) produces no gradient and hence no flux.
void test_conduction_row_sums() {
  const std::vector<Vec3> coords = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                    {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
  const std::vector<Real> Kt = fem::element_conduction(ElementType::C3D8, coords, 50.0);
  for (int a = 0; a < 8; ++a) {
    Real row = 0.0;
    for (int b = 0; b < 8; ++b) row += Kt[static_cast<std::size_t>(a * 8 + b)];
    CX_NEAR(row, 0.0, 1e-9);
  }
  // Symmetry.
  for (int a = 0; a < 8; ++a)
    for (int b = 0; b < 8; ++b)
      CX_NEAR(Kt[static_cast<std::size_t>(a * 8 + b)],
              Kt[static_cast<std::size_t>(b * 8 + a)], 1e-12);
}

// Capacitance of a unit cube sums to rho*c*V (partition of unity: Σ_a Σ_b C[a][b] =
// rho*c ∫ (Σ N_a)(Σ N_b) dV = rho*c*V).
void test_capacitance_total() {
  const std::vector<Vec3> coords = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                    {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
  const std::vector<Real> Ce = fem::element_capacitance(ElementType::C3D8, coords, 3.0);
  const Real total = std::accumulate(Ce.begin(), Ce.end(), 0.0);
  CX_NEAR(total, 3.0, 1e-9);  // rho*c=3, V=1
}

// Analytical 1-D bar: unit cube, x=0 face at T=0, x=1 face at T=100, k=50. The
// steady field is T = 100 x, so every node's temperature is 100*x, and the total
// heat rate through the bar is Q = k A dT / L = 50 * 1 * 100 / 1 = 5000, recovered
// as the sum of the flux reactions on the hot (or cold) face.
void test_bar_linear_profile() {
  Model m = unit_cube_hex(50.0);
  for (Index nd : {1, 4, 5, 8}) m.temp_bcs.push_back(TempBc{nd, 0.0, ""});    // x=0
  for (Index nd : {2, 3, 6, 7}) m.temp_bcs.push_back(TempBc{nd, 100.0, ""});  // x=1

  const ThermalFields t = numerics::solve_heat_transfer(m);
  // Node index i corresponds to node id i+1; check T = 100*x.
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i) {
    const Real x = m.mesh.nodes()[i].x[0];
    CX_NEAR(t.temperature[i], 100.0 * x, 1e-6);
  }
  // Energy balance: reactions sum to zero overall; the cold face (x=0) absorbs -Q,
  // the hot face (x=1) supplies +Q.
  Real cold = 0.0, hot = 0.0, all = 0.0;
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i) {
    all += t.flux_reaction[i];
    if (m.mesh.nodes()[i].x[0] == 0.0) cold += t.flux_reaction[i];
    else hot += t.flux_reaction[i];
  }
  CX_NEAR(all, 0.0, 1e-6);
  CX_NEAR(hot, 5000.0, 1e-6);
  CX_NEAR(cold, -5000.0, 1e-6);
}

// *CFLUX puts the concentrated value straight onto the node's thermal rhs.
void test_cflux_vector() {
  Model m = unit_cube_hex(50.0);
  m.cfluxes.push_back(Cflux{7, 42.0, ""});
  const std::vector<Real> q = fem::thermal_load_vector(m);
  const Index n7 = m.mesh.node_index(7);
  CX_NEAR(q[static_cast<std::size_t>(n7)], 42.0, 1e-12);
  Real total = std::accumulate(q.begin(), q.end(), 0.0);
  CX_NEAR(total, 42.0, 1e-12);  // only one loaded node
}

// *DFLUX over an element face integrates N over the face: the consistent nodal
// fluxes sum to flux * face_area. Face 1 of the unit hex (z=0) has area 1.
void test_dflux_vector() {
  Model m = unit_cube_hex(50.0);
  m.dfluxes.push_back(Dflux{1, 1, 10.0, ""});  // S1, flux 10 -> total 10*Area = 10
  const std::vector<Real> q = fem::thermal_load_vector(m);
  const Real total = std::accumulate(q.begin(), q.end(), 0.0);
  CX_NEAR(total, 10.0, 1e-9);
  // All the flux lands on the four z=0 face nodes (1,2,3,4); interior/top get none.
  for (Index nd : {5, 6, 7, 8}) {
    const Index ni = m.mesh.node_index(nd);
    CX_NEAR(q[static_cast<std::size_t>(ni)], 0.0, 1e-12);
  }
}

// End to end through the parser: a *HEAT TRANSFER, STEADY STATE deck with prescribed
// temperatures on DOF 11 solves the linear bar. Exercises the parser routing of
// *CONDUCTIVITY, *BOUNDARY dof 11, and procedure auto-detection.
void test_parser_and_solve() {
  const std::string deck = R"(
*NODE, NSET=Nall
1, 0., 0., 0.
2, 1., 0., 0.
3, 1., 1., 0.
4, 0., 1., 0.
5, 0., 0., 1.
6, 1., 0., 1.
7, 1., 1., 1.
8, 0., 1., 1.
*ELEMENT, TYPE=C3D8, ELSET=Eall
1, 1, 2, 3, 4, 5, 6, 7, 8
*NSET, NSET=COLD
1, 4, 5, 8
*NSET, NSET=HOT
2, 3, 6, 7
*MATERIAL, NAME=EL
*CONDUCTIVITY
50.
*SOLID SECTION, ELSET=Eall, MATERIAL=EL
*STEP
*HEAT TRANSFER, STEADY STATE
*BOUNDARY
COLD, 11, 11, 0.
HOT, 11, 11, 100.
*END STEP
)";
  const Model m = io::parse_inp(deck);
  CX_CHECK(m.procedure == Procedure::HeatTransferSteady);
  CX_CHECK(m.temp_bcs.size() == 8);
  const ThermalFields t = numerics::solve_heat_transfer(m);
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i)
    CX_NEAR(t.temperature[i], 100.0 * m.mesh.nodes()[i].x[0], 1e-6);
}

// Parser routing for the Phase-3 thermal cards: *FILM / *RADIATE / *PHYSICAL
// CONSTANTS / *INITIAL CONDITIONS, and transient dispatch (a *HEAT TRANSFER step
// without STEADY STATE is now transient, not an error).
void test_parser_film_radiate_transient() {
  const std::string deck = R"(
*NODE, NSET=Nall
1, 0., 0., 0.
2, 1., 0., 0.
3, 1., 1., 0.
4, 0., 1., 0.
5, 0., 0., 1.
6, 1., 0., 1.
7, 1., 1., 1.
8, 0., 1., 1.
*ELEMENT, TYPE=C3D8, ELSET=Eall
1, 1, 2, 3, 4, 5, 6, 7, 8
*MATERIAL, NAME=EL
*CONDUCTIVITY
50.
*DENSITY
7800.
*SPECIFIC HEAT
500.
*PHYSICAL CONSTANTS, ABSOLUTE ZERO=-273.15, STEFAN BOLTZMANN=5.669E-8
*SOLID SECTION, ELSET=Eall, MATERIAL=EL
*INITIAL CONDITIONS, TYPE=TEMPERATURE
Nall, 25.
*STEP
*HEAT TRANSFER
1., 10.
*FILM
1, F4, -12., 10.
*RADIATE
1, R4, 300., 0.8
*END STEP
)";
  const Model m = io::parse_inp(deck);
  CX_CHECK(m.procedure == Procedure::HeatTransferTransient);
  CX_CHECK(m.films.size() == 1);
  CX_CHECK(m.films[0].face == 4);
  CX_NEAR(m.films[0].sink_temp, -12.0, 1e-12);
  CX_NEAR(m.films[0].h, 10.0, 1e-12);
  CX_CHECK(m.radiates.size() == 1);
  CX_CHECK(m.radiates[0].face == 4);
  CX_NEAR(m.radiates[0].emissivity, 0.8, 1e-12);
  CX_NEAR(m.physical.absolute_zero, -273.15, 1e-9);
  CX_NEAR(m.physical.sigma, 5.669e-8, 1e-15);
  CX_NEAR(m.initial_temperature, 25.0, 1e-12);
  CX_NEAR(m.increment.total, 10.0, 1e-12);
  CX_NEAR(m.increment.initial, 1.0, 1e-12);
}

// Integration-point heat flux HFL (task 6.2). The linear bar (T = 100 x, k = 50) has
// a UNIFORM temperature gradient dT/dx = 100, so q = -k grad(T) gives qx = -5000 at
// EVERY Gauss point (transverse components zero), exactly the Fourier flux q = k dT/L
// (here L = 1). We check both the raw integration-point HFL and the extrapolated
// nodal HFL against that exact value.
void test_heat_flux_bar() {
  Model m = unit_cube_hex(50.0);
  for (Index nd : {1, 4, 5, 8}) m.temp_bcs.push_back(TempBc{nd, 0.0, ""});    // x=0
  for (Index nd : {2, 3, 6, 7}) m.temp_bcs.push_back(TempBc{nd, 100.0, ""});  // x=1

  const ThermalFields t = numerics::solve_heat_transfer(m);
  CX_CHECK(!t.hfl_points.empty());
  for (const HeatFluxPoint& p : t.hfl_points) {
    CX_CHECK(p.elem_id == 1);
    CX_NEAR(p.flux[0], -5000.0, 1e-6);  // qx = -k dT/dx = -50*100
    CX_NEAR(p.flux[1], 0.0, 1e-6);
    CX_NEAR(p.flux[2], 0.0, 1e-6);
  }
  // Nodal (extrapolated) HFL is the same uniform vector at every node.
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i) {
    CX_NEAR(t.heat_flux[i][0], -5000.0, 1e-6);
    CX_NEAR(t.heat_flux[i][1], 0.0, 1e-6);
    CX_NEAR(t.heat_flux[i][2], 0.0, 1e-6);
  }
}

// PropertyTable interpolation: piecewise-linear with clamping outside the range, and
// a single-row table is a constant.
void test_property_table_interpolation() {
  PropertyTable single{7.5};
  CX_NEAR(single.at(-100.0), 7.5, 1e-12);
  CX_NEAR(single.at(1000.0), 7.5, 1e-12);

  PropertyTable table;
  table.value = {10.0, 20.0, 50.0};
  table.temp = {0.0, 100.0, 200.0};
  CX_NEAR(table.at(-50.0), 10.0, 1e-12);   // clamp below
  CX_NEAR(table.at(0.0), 10.0, 1e-12);
  CX_NEAR(table.at(50.0), 15.0, 1e-12);    // midway 10..20
  CX_NEAR(table.at(100.0), 20.0, 1e-12);
  CX_NEAR(table.at(150.0), 35.0, 1e-12);   // midway 20..50
  CX_NEAR(table.at(200.0), 50.0, 1e-12);
  CX_NEAR(table.at(500.0), 50.0, 1e-12);   // clamp above
}

// Temperature-dependent conductivity k(T) = k0 + s (T - Tref) on the 1-D bar. For a
// 1-D steady bar with no source the flux q is CONSTANT along x, so the Kirchhoff
// transform gives  q L = ∫_{T0}^{T1} k(T) dT = k0 (T1-T0) + (s/2)((T1-Tref)^2 -
// (T0-Tref)^2).  With T0=0, T1=100, k0=50, s=0.5, Tref=0:
//   ∫ = 50*100 + 0.25*(100^2 - 0) = 5000 + 2500 = 7500, so |q| = 7500 (L=1, A=1).
// The temperature PROFILE is nonlinear (the mid node is NOT at 50). We check the FEM
// (Picard-iterated) matches: the total heat rate through each face and the
// Kirchhoff-implied mid-plane temperature.
void test_temp_dependent_conductivity_bar() {
  // A two-hex bar along x (0..2) so there is an interior mid-plane node to check the
  // nonlinear profile; unit cross section. Nodes: plane x=0 {1,4,5,8}, x=1
  // {2,3,6,7}, x=2 {9,10,11,12} — build explicitly.
  Model m;
  m.mesh.add_node(1, {0, 0, 0});
  m.mesh.add_node(2, {1, 0, 0});
  m.mesh.add_node(3, {1, 1, 0});
  m.mesh.add_node(4, {0, 1, 0});
  m.mesh.add_node(5, {0, 0, 1});
  m.mesh.add_node(6, {1, 0, 1});
  m.mesh.add_node(7, {1, 1, 1});
  m.mesh.add_node(8, {0, 1, 1});
  m.mesh.add_node(9, {2, 0, 0});
  m.mesh.add_node(10, {2, 1, 0});
  m.mesh.add_node(11, {2, 0, 1});
  m.mesh.add_node(12, {2, 1, 1});
  m.mesh.add_element(1, ElementType::C3D8, {1, 2, 3, 4, 5, 6, 7, 8});
  m.mesh.add_element(2, ElementType::C3D8, {2, 9, 10, 3, 6, 11, 12, 7});
  m.mesh.add_elset("EALL", {1, 2});
  Material mat;
  mat.name = "EL";
  Thermal th;
  th.conductivity.value = {50.0, 100.0};  // k0=50 at T=0, 100 at T=100 -> s=0.5
  th.conductivity.temp = {0.0, 100.0};
  mat.thermal = th;
  m.materials["EL"] = mat;
  m.sections.push_back(SolidSection{"EALL", "EL"});
  m.procedure = Procedure::HeatTransferSteady;
  for (Index nd : {1, 4, 5, 8}) m.temp_bcs.push_back(TempBc{nd, 0.0, ""});      // x=0
  for (Index nd : {9, 10, 11, 12}) m.temp_bcs.push_back(TempBc{nd, 100.0, ""}); // x=2

  const ThermalFields t = numerics::solve_heat_transfer(m);

  // Kirchhoff: q L = ∫_0^100 (50 + 0.5 T) dT = 5000 + 2500 = 7500 with L=2, A=1, so
  // the constant flux is q = 7500 / L = 3750 and the heat rate Q = q A = 3750, the
  // sum of the flux reactions on the hot (x=2) face.
  Real hot = 0.0;
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i)
    if (m.mesh.nodes()[i].x[0] == 2.0) hot += t.flux_reaction[i];
  CX_NEAR(hot, 3750.0, 1e-2);

  // Mid-plane (x=1) temperature Tm solves ∫_0^{Tm} k dT = (x/L) ∫_0^{T1} k dT with
  // x/L = 1/2: 50 Tm + 0.25 Tm^2 = 3750 -> Tm^2 + 200 Tm - 15000 = 0 ->
  // Tm = -100 + sqrt(10000 + 15000) = -100 + sqrt(25000) = 58.113883...
  const Real Tm = -100.0 + std::sqrt(25000.0);
  const Index n2 = m.mesh.node_index(2);
  CX_NEAR(t.temperature[static_cast<std::size_t>(n2)], Tm, 1e-3);
  // A CONSTANT conductivity would put the mid node at 50; confirm we are meaningfully
  // above that (nonlinear k(T) really took effect).
  CX_CHECK(t.temperature[static_cast<std::size_t>(n2)] > 55.0);
}

// A single-row *CONDUCTIVITY table is byte-for-byte the constant path: k(T) resolves
// to the same conductivity regardless of temperature, so the linear-bar solve is
// unchanged and has_temp_dependent_thermal() is false.
void test_constant_table_unchanged() {
  Model m = unit_cube_hex(50.0);  // single-row Thermal{50, 0}
  CX_CHECK(!m.has_temp_dependent_thermal());
  for (Index nd : {1, 4, 5, 8}) m.temp_bcs.push_back(TempBc{nd, 0.0, ""});
  for (Index nd : {2, 3, 6, 7}) m.temp_bcs.push_back(TempBc{nd, 100.0, ""});
  const ThermalFields t = numerics::solve_heat_transfer(m);
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i)
    CX_NEAR(t.temperature[i], 100.0 * m.mesh.nodes()[i].x[0], 1e-9);
}

// Parser: a multi-row *CONDUCTIVITY / *SPECIFIC HEAT / *EXPANSION builds a table; a
// non-increasing temperature column is rejected.
void test_parser_temp_tables() {
  const std::string deck = R"(
*MATERIAL, NAME=EL
*CONDUCTIVITY
50., 0.
100., 200.
*SPECIFIC HEAT
400., 0.
600., 100.
*EXPANSION, ZERO=20.
1.0e-5, 0.
1.4e-5, 500.
)";
  const Model m = io::parse_inp(deck);
  const Material& mat = m.materials.at("EL");
  CX_CHECK(mat.thermal->conductivity.value.size() == 2);
  CX_NEAR(mat.thermal->conductivity.at(100.0), 75.0, 1e-9);   // midway 50..100
  CX_CHECK(mat.thermal->specific_heat.value.size() == 2);
  CX_NEAR(mat.thermal->specific_heat.at(50.0), 500.0, 1e-9);
  CX_CHECK(mat.expansion->alpha.value.size() == 2);
  CX_NEAR(mat.expansion->t_ref, 20.0, 1e-12);
  CX_NEAR(mat.expansion->alpha.at(250.0), 1.2e-5, 1e-12);
  CX_CHECK(m.has_temp_dependent_thermal());

  bool threw = false;
  try {
    io::parse_inp("*MATERIAL, NAME=X\n*CONDUCTIVITY\n50., 100.\n60., 100.\n");
  } catch (const io::ParseError&) {
    threw = true;
  }
  CX_CHECK(threw);
}

}  // namespace

int main() {
  test_conduction_row_sums();
  test_capacitance_total();
  test_bar_linear_profile();
  test_cflux_vector();
  test_dflux_vector();
  test_parser_and_solve();
  test_parser_film_radiate_transient();
  test_heat_flux_bar();
  test_property_table_interpolation();
  test_temp_dependent_conductivity_bar();
  test_constant_table_unchanged();
  test_parser_temp_tables();
  if (cxtest::g_failures == 0) std::printf("test_thermal: OK\n");
  CX_MAIN_RETURN();
}
