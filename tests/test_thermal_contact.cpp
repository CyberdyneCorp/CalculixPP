// Thermal contact — *GAP CONDUCTANCE / *GAP HEAT GENERATION (spec: contact — thermal
// contact conductance and gap heat generation; tasks.md 2.4). Validation is ANALYTICAL:
// two unit cubes stacked in z that touch ONLY through a contact interface carrying a gap
// conductance, so all heat between them crosses the gap by q = h_c (T_slave - T_master).
//
// SERIES-RESISTANCE MODEL: cube 1 spans z in [0,1] (its base z=0 held at T_cold), cube 2
// z in [1,2] (its top z=2 held at T_hot), with distinct nodes on z=1 coupled only by the
// *GAP CONDUCTANCE h. Each unit cube of conductivity k, unit cross-section, unit length,
// is a thermal resistor R_cube = L/(k A) = 1/k; the gap adds R_gap = 1/(h A) = 1/h. The
// three resistors are in series, so the through-flux is
//   Q = (T_hot - T_cold) / (R_cube + R_gap + R_cube) = dT / (2/k + 1/h),
// and the temperature DROP across the interface is
//   dT_gap = Q * R_gap = Q / h,
// i.e. the conductance-controlled interface flux q = h dT_gap = Q. This is the defining
// gap-conductance relation. We check Q (energy through the stack) and the interface drop.
#include <cmath>
#include <string>
#include <vector>

#include "calculixpp/io/inp_parser.hpp"
#include "calculixpp/numerics/heat_transfer.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

// Two unit C3D8 cubes stacked in z, distinct nodes on the z=1 interface. Cube 1 nodes
// 1-8 (z=0 then z=1), cube 2 nodes 9-16 (z=1 then z=2). Conductivity `k`, gap conductance
// `h`. Bottom face (z=0, nodes 1-4) held at `t_cold`; top face (z=2, nodes 13-16) held at
// `t_hot`. Slave = top cube bottom nodes (9-12) as a node surface; master = bottom cube
// top face S2. Steady heat transfer.
std::string gap_deck(double k, double h, double t_cold, double t_hot) {
  return std::string(R"(
*NODE, NSET=NALL
1, 0., 0., 0.
2, 1., 0., 0.
3, 1., 1., 0.
4, 0., 1., 0.
5, 0., 0., 1.
6, 1., 0., 1.
7, 1., 1., 1.
8, 0., 1., 1.
9, 0., 0., 1.
10, 1., 0., 1.
11, 1., 1., 1.
12, 0., 1., 1.
13, 0., 0., 2.
14, 1., 0., 2.
15, 1., 1., 2.
16, 0., 1., 2.
*ELEMENT, TYPE=C3D8, ELSET=EBASE
1, 1,2,3,4,5,6,7,8
*ELEMENT, TYPE=C3D8, ELSET=ETOP
2, 9,10,11,12,13,14,15,16
*ELSET, ELSET=EALL
EBASE, ETOP
*NSET, NSET=NCOLD
1,2,3,4
*NSET, NSET=NHOT
13,14,15,16
*MATERIAL, NAME=STEEL
*ELASTIC
210000., 0.3
*CONDUCTIVITY
)") + std::to_string(k) + R"(
*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL
*SURFACE, NAME=SMAST
EBASE, S2
*SURFACE, NAME=SSLAV, TYPE=NODE
9,10,11,12
*SURFACE INTERACTION, NAME=SI
*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=LINEAR
1.E7
*GAP CONDUCTANCE
)" + std::to_string(h) + R"(
*CONTACT PAIR, INTERACTION=SI, TYPE=NODE TO SURFACE
SSLAV, SMAST
*STEP
*HEAT TRANSFER, STEADY STATE
*BOUNDARY
NCOLD,11,11,)" + std::to_string(t_cold) + R"(
NHOT,11,11,)" + std::to_string(t_hot) + R"(
*NODE PRINT, NSET=NALL
NT
*END STEP
)";
}

double node_T(const ThermalFields& t, const Model& m, Index id) {
  return t.temperature[static_cast<std::size_t>(m.mesh.node_index(id))];
}

// The through-flux from the RFL reaction at the cold face (sum of the heat-flux reaction
// over the held-cold nodes) equals the series flux Q, and the interface temperature drop
// equals Q/h.
void test_gap_conductance_series_resistance() {
  const double k = 1.0, h = 1.0, tc = 0.0, thot = 100.0;
  const Model m = io::parse_inp(gap_deck(k, h, tc, thot));
  CX_CHECK(m.has_thermal_contact());
  const ThermalFields t = numerics::solve_heat_transfer(m);

  // Analytical series flux Q = dT / (2/k + 1/h). A = 1.
  const double Q_expected = (thot - tc) / (2.0 / k + 1.0 / h);

  // Interface temperatures: master (bottom-cube top, nodes 5-8) and slave (top-cube
  // bottom, nodes 9-12). Each face is isothermal by symmetry.
  double tm = 0.0, ts = 0.0;
  for (Index id : {5, 6, 7, 8}) tm += node_T(t, m, id);
  for (Index id : {9, 10, 11, 12}) ts += node_T(t, m, id);
  tm /= 4.0;
  ts /= 4.0;
  const double dT_gap = ts - tm;  // hot side (slave, top) is warmer than master (bottom)

  // Interface flux q = h * A * dT_gap must equal the series flux Q.
  const double q_interface = h * 1.0 * dT_gap;
  CX_NEAR(q_interface, Q_expected, 1e-6 * Q_expected);

  // The interface drop matches Q/h (the gap-resistance temperature drop).
  CX_NEAR(dT_gap, Q_expected / h, 1e-6 * (Q_expected / h));

  // The cold-face RFL reaction (heat entering the stack) balances the through-flux Q.
  double rfl = 0.0;
  for (Index id : {1, 2, 3, 4})
    rfl += t.flux_reaction[static_cast<std::size_t>(m.mesh.node_index(id))];
  CX_NEAR(std::fabs(rfl), Q_expected, 1e-6 * Q_expected);
}

// A stiffer gap conductance (larger h) shrinks the interface temperature drop toward zero
// (perfect thermal contact); a softer one widens it. In the limit h -> inf the interface
// drop -> 0 and Q -> the pure-conduction flux dT/(2/k).
void test_conductance_controls_interface_drop() {
  const double k = 1.0, tc = 0.0, thot = 100.0;
  auto interface_drop = [&](double h) {
    const Model m = io::parse_inp(gap_deck(k, h, tc, thot));
    const ThermalFields t = numerics::solve_heat_transfer(m);
    double tm = 0.0, ts = 0.0;
    for (Index id : {5, 6, 7, 8}) tm += node_T(t, m, id);
    for (Index id : {9, 10, 11, 12}) ts += node_T(t, m, id);
    return (ts - tm) / 4.0;
  };
  const double soft = interface_drop(0.5);   // large gap resistance -> big drop
  const double stiff = interface_drop(100.0); // near-perfect contact -> tiny drop
  CX_CHECK(soft > stiff);
  CX_CHECK(stiff < 2.0);            // near-perfect contact almost eliminates the drop
  CX_CHECK(soft > 40.0);            // R_gap = 2 dominates R_cube = 1 -> ~half of dT
}

// Gap heat generation: with both faces held at the SAME temperature (no conduction drop)
// a *GAP HEAT GENERATION q deposits heat q*A at the interface, split evenly. With the
// stack's two ends held equal, the deposited heat must flow out both ends: the total
// RFL over both held faces balances q*A (energy conservation), confirming the source.
void test_gap_heat_generation_energy_balance() {
  const double k = 1.0, h = 1.0, q_gen = 20.0, tset = 50.0;
  std::string deck = std::string(R"(
*NODE, NSET=NALL
1, 0., 0., 0.
2, 1., 0., 0.
3, 1., 1., 0.
4, 0., 1., 0.
5, 0., 0., 1.
6, 1., 0., 1.
7, 1., 1., 1.
8, 0., 1., 1.
9, 0., 0., 1.
10, 1., 0., 1.
11, 1., 1., 1.
12, 0., 1., 1.
13, 0., 0., 2.
14, 1., 0., 2.
15, 1., 1., 2.
16, 0., 1., 2.
*ELEMENT, TYPE=C3D8, ELSET=EBASE
1, 1,2,3,4,5,6,7,8
*ELEMENT, TYPE=C3D8, ELSET=ETOP
2, 9,10,11,12,13,14,15,16
*ELSET, ELSET=EALL
EBASE, ETOP
*NSET, NSET=NCOLD
1,2,3,4
*NSET, NSET=NHOT
13,14,15,16
*MATERIAL, NAME=STEEL
*ELASTIC
210000., 0.3
*CONDUCTIVITY
)") + std::to_string(k) + R"(
*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL
*SURFACE, NAME=SMAST
EBASE, S2
*SURFACE, NAME=SSLAV, TYPE=NODE
9,10,11,12
*SURFACE INTERACTION, NAME=SI
*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=LINEAR
1.E7
*GAP CONDUCTANCE
)" + std::to_string(h) + R"(
*GAP HEAT GENERATION
)" + std::to_string(q_gen) + R"(
*CONTACT PAIR, INTERACTION=SI, TYPE=NODE TO SURFACE
SSLAV, SMAST
*STEP
*HEAT TRANSFER, STEADY STATE
*BOUNDARY
NCOLD,11,11,)" + std::to_string(tset) + R"(
NHOT,11,11,)" + std::to_string(tset) + R"(
*NODE PRINT, NSET=NALL
NT
*END STEP
)";
  const Model m = io::parse_inp(deck);
  CX_CHECK(m.has_thermal_contact());
  const ThermalFields t = numerics::solve_heat_transfer(m);

  // Total heat leaving through both held faces balances the generated q*A (A = 1).
  double rfl = 0.0;
  for (Index id : {1, 2, 3, 4, 13, 14, 15, 16})
    rfl += t.flux_reaction[static_cast<std::size_t>(m.mesh.node_index(id))];
  // The generated heat leaves the domain: reactions sum to -q*A (heat removed at the BCs).
  CX_NEAR(std::fabs(rfl), q_gen * 1.0, 1e-6 * q_gen);
}

}  // namespace

int main() {
  test_gap_conductance_series_resistance();
  test_conductance_controls_interface_drop();
  test_gap_heat_generation_energy_balance();
  if (cxtest::g_failures == 0) std::printf("test_thermal_contact: OK\n");
  CX_MAIN_RETURN();
}
