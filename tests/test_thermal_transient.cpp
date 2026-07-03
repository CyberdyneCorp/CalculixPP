// Transient conduction, convective film, and surface radiation
// (spec: heat-transfer-analysis, tasks 3.2 / 3.3 / 3.4).
//
// Validation strategy (self-contained analytical cases; no reference deck needed
// for transient/radiation — the film STEADY case is separately matched against the
// stock CalculiX oneel20fi deck by the Python regression suite):
//
//  * Transient uniform-field lumped cooling: a unit cube with the SAME film on all
//    six faces to a common sink keeps the temperature field spatially UNIFORM (a
//    uniform field is an exact eigenmode — conduction row sums vanish, and the
//    consistent film/capacitance matrices preserve uniformity). The FEM then reduces
//    EXACTLY to the lumped ODE  Cv dT/dt = -h A (T - T_sink), whose backward-Euler
//    recurrence  T_{n+1} = (Cv/dt T_n + h A T_sink)/(Cv/dt + h A)  the driver must
//    reproduce to machine precision, and which tracks the exact exponential
//    T(t) = T_sink + (T0 - T_sink) e^{-hA/Cv t} within the time-discretization error.
//
//  * Transient -> steady limit: with a long period the transient field converges to
//    the independently-computed steady conduction+film solution.
//
//  * Surface radiation equilibrium: a bar with a fixed hot face and a radiating cold
//    face reaches the temperature where conduction in equals radiation out,
//    k A (T_hot - T)/L = eps sigma A (T_abs^4 - T_amb_abs^4), checked against the
//    scalar root of that balance.
#include <cmath>
#include <vector>

#include "calculixpp/fem/element.hpp"
#include "calculixpp/numerics/heat_transfer.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

// Unit cube as one C3D8 hex, CalculiX node order, isotropic conductivity k.
Model unit_cube(Real k) {
  Model m;
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
  return m;
}

// Uniform-field lumped cooling. The cube starts at T0 and all six faces convect to a
// sink; the field stays uniform so the FEM equals the lumped ODE exactly.
void test_transient_uniform_lumped() {
  const Real k = 1.0, rho = 2.0, c = 3.0, h = 5.0, T0 = 100.0, T_sink = 20.0;
  const Real Cv = rho * c;         // volume = 1 -> total capacitance rho*c
  const Real A = 6.0;              // six unit faces
  Model m = unit_cube(k);
  m.materials["EL"].density = rho;
  m.materials["EL"].thermal->specific_heat = c;
  m.procedure = Procedure::HeatTransferTransient;
  for (int face = 1; face <= 6; ++face)
    m.films.push_back(Film{1, face, T_sink, h, ""});
  m.initial_temperature = T0;

  const Real dt = 0.01, period = 0.2;
  m.increment.initial = dt;
  m.increment.total = period;

  const ThermalFields t = numerics::solve_heat_transfer(m);

  // Exact backward-Euler recurrence for the lumped mode, replayed here. The film
  // sink carries no *AMPLITUDE, so (like every unamplified load) it ramps linearly
  // 0 -> T_sink over the step: at step s+1 the sink fraction is lambda = (s+1)/steps.
  const int steps = static_cast<int>(std::llround(period / dt));
  Real T = T0;
  for (int s = 0; s < steps; ++s) {
    const Real lambda = static_cast<Real>(s + 1) / static_cast<Real>(steps);
    const Real sink = lambda * T_sink;
    T = (Cv / dt * T + h * A * sink) / (Cv / dt + h * A);
  }

  // Every node holds the same (uniform) temperature = the lumped recurrence value.
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i)
    CX_NEAR(t.temperature[i], T, 1e-8);

  CX_CHECK(T < T0 && T > T_sink);        // cooled toward the (ramping) sink, monotone
}

// A long transient relaxes to the steady conduction+film equilibrium. Cube with the
// x=0 face fixed at 0 and a film on the x=1 face (face S3 in hex numbering carries
// the +x face) reaches the steady solve's field.
void test_transient_to_steady() {
  const Real k = 50.0, rho = 1.0, c = 1.0, h = 10.0;
  Model m = unit_cube(k);
  m.materials["EL"].density = rho;
  m.materials["EL"].thermal->specific_heat = c;
  for (Index nd : {1, 4, 5, 8}) m.temp_bcs.push_back(TempBc{nd, 0.0, ""});  // x=0 face
  // +x face is hex face 4 (nodes 2,6,7,3) -> convect to sink -12, h=10.
  m.films.push_back(Film{1, 4, -12.0, h, ""});

  // Steady reference.
  Model ms = m;
  ms.procedure = Procedure::HeatTransferSteady;
  const ThermalFields steady = numerics::solve_heat_transfer(ms);

  // Transient over a long period from a zero initial field.
  m.procedure = Procedure::HeatTransferTransient;
  m.initial_temperature = 0.0;
  m.increment.initial = 0.05;
  m.increment.total = 50.0;  // long enough to relax (sink fully ramped, field settled)
  const ThermalFields tr = numerics::solve_heat_transfer(m);

  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i)
    CX_NEAR(tr.temperature[i], steady.temperature[i], 1e-3);
}

// Surface radiation equilibrium: hot face fixed at T_hot, opposite face radiates to
// ambient. The steady field is uniform along x per node group; the +x face reaches
// the balance k(T_hot - T)/L = eps sigma (T_abs^4 - T_amb_abs^4) with L=1, A=1.
void test_radiation_equilibrium() {
  const Real k = 50.0, eps = 0.8, sigma = 5.669e-8, az = -273.15;
  const Real T_hot = 1000.0, T_amb = 300.0;
  Model m = unit_cube(k);
  m.procedure = Procedure::HeatTransferSteady;
  m.physical = PhysicalConstants{az, sigma};
  for (Index nd : {1, 4, 5, 8}) m.temp_bcs.push_back(TempBc{nd, T_hot, ""});  // x=0 hot
  m.radiates.push_back(Radiate{1, 4, T_amb, eps});  // +x face radiates

  const ThermalFields t = numerics::solve_heat_transfer(m);

  // Cold-face temperature (node 2, x=1). Solve the 1-D balance for the reference:
  //   k (T_hot - T)/L = eps sigma ((T-az)^4 - (T_amb-az)^4),  L = 1.
  auto balance = [&](Real T) {
    const Real Ta = T - az, Tamb = T_amb - az;
    return k * (T_hot - T) - eps * sigma * (Ta * Ta * Ta * Ta - Tamb * Tamb * Tamb * Tamb);
  };
  Real lo = T_amb, hi = T_hot;  // monotone-decreasing residual in T
  for (int i = 0; i < 200; ++i) {
    const Real mid = 0.5 * (lo + hi);
    if (balance(mid) > 0.0) lo = mid; else hi = mid;
  }
  const Real T_ref = 0.5 * (lo + hi);

  const Index n2 = m.mesh.node_index(2);  // x=1 face corner
  CX_NEAR(t.temperature[static_cast<std::size_t>(n2)], T_ref, 1e-3);
}

}  // namespace

int main() {
  test_transient_uniform_lumped();
  test_transient_to_steady();
  test_radiation_equilibrium();
  if (cxtest::g_failures == 0) std::printf("test_thermal_transient: OK\n");
  CX_MAIN_RETURN();
}
