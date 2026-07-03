// Thermal-expansion coupling and one-way coupled temperature-displacement
// (spec: heat-transfer-analysis — coupled; tasks 4.1 / 4.3).
//
// Validated ANALYTICALLY (no reference deck needed):
//   (a) Free bar: a bar heated uniformly by dT with no axial constraint expands by
//       alpha*dT*L and carries ZERO stress (eps_mech = eps_th everywhere).
//   (b) Axially-constrained bar: a bar with u_x = 0 on both end faces (lateral faces
//       free) heated by dT develops the uniaxial thermal stress sigma_xx = -E*alpha*dT
//       (and sigma_yy = sigma_zz = 0, since the lateral strains relax freely).
//   (c) One-way *COUPLED TEMPERATURE-DISPLACEMENT end to end: a bar with a prescribed
//       temperature difference is solved thermally (linear T profile), the resulting
//       thermal strain is applied, and the axially-constrained mechanical solve gives
//       sigma_xx = -E*alpha*T_bar with T_bar the mean temperature rise.
#include <cmath>
#include <string>
#include <vector>

#include "calculixpp/fem/loads.hpp"
#include "calculixpp/io/inp_parser.hpp"
#include "calculixpp/numerics/heat_transfer.hpp"
#include "calculixpp/numerics/linear_static.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

constexpr Real kE = 210000.0;
constexpr Real kNu = 0.3;
constexpr Real kAlpha = 1.2e-5;

// A unit cube meshed as one C3D8 hex, x the bar axis. Linear-elastic with *EXPANSION.
Model unit_cube_mech() {
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
  mat.elastic = ElasticIso{kE, kNu};
  mat.expansion = Expansion{kAlpha, 0.0};  // Tref = 0
  m.materials["EL"] = mat;
  m.sections.push_back(SolidSection{"EALL", "EL"});
  m.procedure = Procedure::Static;
  return m;
}

Real von_mises(const Voigt6& s) {
  const Real a = s[0] - s[1], b = s[1] - s[2], c = s[2] - s[0];
  return std::sqrt(0.5 * (a * a + b * b + c * c) +
                   3.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]));
}

// (a) Free bar: uniform heating dT, only rigid-body motion removed. Expands by
// alpha*dT along every axis (isotropic), zero stress.
void test_free_bar_expands_stress_free() {
  Model m = unit_cube_mech();
  const Real dT = 100.0;
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i)
    m.applied_temperature[m.mesh.nodes()[i].id] = dT;
  // Statically-determinate rigid-body suppression: pin node 1 fully, node 2 in y,z,
  // node 4 in z — leaves all faces free to expand isotropically.
  m.spcs.push_back(Spc{1, 1, 0.0, ""});
  m.spcs.push_back(Spc{1, 2, 0.0, ""});
  m.spcs.push_back(Spc{1, 3, 0.0, ""});
  m.spcs.push_back(Spc{2, 2, 0.0, ""});
  m.spcs.push_back(Spc{2, 3, 0.0, ""});
  m.spcs.push_back(Spc{4, 3, 0.0, ""});

  CX_CHECK(m.has_thermal_strain());
  const StaticFields f = numerics::solve_linear_static(m);

  // Node 7 (1,1,1) should displace by alpha*dT along each axis (free thermal growth).
  const Index n7 = m.mesh.node_index(7);
  const Real expect = kAlpha * dT;  // *L=1
  for (int c = 0; c < 3; ++c)
    CX_NEAR(f.displacement[static_cast<std::size_t>(n7)][static_cast<std::size_t>(c)],
            expect, 1e-9);
  // Stress is (essentially) zero everywhere: eps_mech = eps_th.
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i)
    CX_NEAR(von_mises(f.stress[i]), 0.0, 1e-4);
  // Reactions vanish (no external load, thermally self-equilibrated).
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i)
    for (int c = 0; c < 3; ++c)
      CX_NEAR(f.reaction[i][static_cast<std::size_t>(c)], 0.0, 1e-4);
}

// (b) Axially-constrained bar: u_x = 0 on both end faces (x=0 and x=1); lateral faces
// free. Uniform heating dT -> uniaxial thermal stress sigma_xx = -E*alpha*dT.
void test_constrained_bar_thermal_stress() {
  Model m = unit_cube_mech();
  const Real dT = 100.0;
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i)
    m.applied_temperature[m.mesh.nodes()[i].id] = dT;
  // Fix u_x on both end faces (x=0: nodes 1,4,5,8; x=1: nodes 2,3,6,7).
  for (Index nd : {1, 4, 5, 8, 2, 3, 6, 7}) m.spcs.push_back(Spc{nd, 1, 0.0, ""});
  // Remove y/z rigid-body: pin node 1 in y,z and node 2 in z.
  m.spcs.push_back(Spc{1, 2, 0.0, ""});
  m.spcs.push_back(Spc{1, 3, 0.0, ""});
  m.spcs.push_back(Spc{2, 3, 0.0, ""});

  const StaticFields f = numerics::solve_linear_static(m);
  const Real sxx_expect = -kE * kAlpha * dT;  // = -252.0
  // Every node: sigma_xx = -E*alpha*dT; sigma_yy = sigma_zz = 0 (lateral relaxation).
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i) {
    CX_NEAR(f.stress[i][0], sxx_expect, 1e-3);
    CX_NEAR(f.stress[i][1], 0.0, 1e-3);
    CX_NEAR(f.stress[i][2], 0.0, 1e-3);
  }
  // Total x-reaction on the x=1 face balances the thermal stress over the face area 1.
  Real rx_hi = 0.0;
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i)
    if (m.mesh.nodes()[i].x[0] == 1.0) rx_hi += f.reaction[i][0];
  CX_NEAR(rx_hi, sxx_expect, 1e-3);  // pushing outward is resisted with -E*alpha*dT*A
}

// (b') Temperature-dependent expansion alpha(T): an axially-constrained bar heated
// UNIFORMLY by dT (so the element mean temperature is dT) develops sigma_xx = -E *
// alpha(dT) * dT with alpha evaluated from the *EXPANSION table at the mean temperature.
// Table: alpha = 1.0e-5 at T=0, 2.0e-5 at T=200, Tref=0; at dT=100 -> alpha=1.5e-5.
void test_constrained_bar_temp_dependent_alpha() {
  Model m = unit_cube_mech();
  Expansion ex;
  ex.alpha.value = {1.0e-5, 2.0e-5};
  ex.alpha.temp = {0.0, 200.0};
  ex.t_ref = 0.0;
  m.materials["EL"].expansion = ex;
  const Real dT = 100.0;
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i)
    m.applied_temperature[m.mesh.nodes()[i].id] = dT;
  for (Index nd : {1, 4, 5, 8, 2, 3, 6, 7}) m.spcs.push_back(Spc{nd, 1, 0.0, ""});
  m.spcs.push_back(Spc{1, 2, 0.0, ""});
  m.spcs.push_back(Spc{1, 3, 0.0, ""});
  m.spcs.push_back(Spc{2, 3, 0.0, ""});

  const StaticFields f = numerics::solve_linear_static(m);
  const Real alpha_at_dT = 1.5e-5;  // interpolated at T=100 in [0,200]
  const Real sxx_expect = -kE * alpha_at_dT * dT;  // = -315.0
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i) {
    CX_NEAR(f.stress[i][0], sxx_expect, 1e-3);
    CX_NEAR(f.stress[i][1], 0.0, 1e-3);
    CX_NEAR(f.stress[i][2], 0.0, 1e-3);
  }
}

// Sanity: with NO applied temperature the thermal path is inert (mechanical unchanged).
void test_no_temperature_is_pure_mechanics() {
  Model m = unit_cube_mech();
  CX_CHECK(!m.has_thermal_strain());  // no applied_temperature
  const std::vector<Real> f_th = fem::thermal_strain_load_vector(m);
  for (const Real v : f_th) CX_NEAR(v, 0.0, 0.0);
}

// (c) One-way *COUPLED TEMPERATURE-DISPLACEMENT end to end via the deck path. A bar
// with hot/cold end faces solves a linear temperature profile; the axially-restrained
// mechanical field then carries sigma_xx = -E*alpha*T_bar, T_bar the mean rise.
void test_coupled_deck_end_to_end() {
  const std::string deck = R"(
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
1, 1, 2, 3, 4, 5, 6, 7, 8
*MATERIAL, NAME=STEEL
*ELASTIC
210000., 0.3
*CONDUCTIVITY
50.
*EXPANSION
1.2e-5
*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL
*STEP
*COUPLED TEMPERATURE-DISPLACEMENT
*BOUNDARY
1, 11, 11, 0.
4, 11, 11, 0.
5, 11, 11, 0.
8, 11, 11, 0.
2, 11, 11, 200.
3, 11, 11, 200.
6, 11, 11, 200.
7, 11, 11, 200.
*BOUNDARY
1, 1
4, 1
5, 1
8, 1
2, 1
3, 1
6, 1
7, 1
1, 2
1, 3
2, 3
*END STEP
)";
  const Model m = io::parse_inp(deck);
  CX_CHECK(m.procedure == Procedure::Coupled);
  // *TEMPERATURE routing: in a coupled/heat step the DOF-11 BCs are thermal, not
  // mechanical SPCs; the mechanical SPCs are only the DOF 1/2/3 boundary lines.
  CX_CHECK(m.temp_bcs.size() == 8);

  const CoupledFields cf = numerics::solve_coupled(m);
  // Thermal: linear profile T = 200*x.
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i)
    CX_NEAR(cf.thermal.temperature[i], 200.0 * m.mesh.nodes()[i].x[0], 1e-6);
  // Mechanical: mean temperature rise T_bar = 100 (Tref=0), so the axially-restrained
  // bar carries the volume-averaged uniaxial thermal stress sigma_xx = -E*alpha*T_bar.
  const Real sxx_expect = -210000.0 * 1.2e-5 * 100.0;  // = -252.0
  // The integrated x-reaction on the hot face balances the mean thermal stress.
  Real rx_hi = 0.0;
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i)
    if (m.mesh.nodes()[i].x[0] == 1.0) rx_hi += cf.mechanical.reaction[i][0];
  CX_NEAR(rx_hi, sxx_expect, 1e-2);
}

// The coupled deck of test_coupled_deck_end_to_end as a Model, for the scheme tests.
Model coupled_bar_model() {
  const std::string deck = R"(
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
1, 1, 2, 3, 4, 5, 6, 7, 8
*MATERIAL, NAME=STEEL
*ELASTIC
210000., 0.3
*CONDUCTIVITY
50.
*EXPANSION
1.2e-5
*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL
*STEP
*COUPLED TEMPERATURE-DISPLACEMENT
*BOUNDARY
1, 11, 11, 0.
4, 11, 11, 0.
5, 11, 11, 0.
8, 11, 11, 0.
2, 11, 11, 200.
3, 11, 11, 200.
6, 11, 11, 200.
7, 11, 11, 200.
*BOUNDARY
1, 1
4, 1
5, 1
8, 1
2, 1
3, 1
6, 1
7, 1
1, 2
1, 3
2, 3
*END STEP
)";
  return io::parse_inp(deck);
}

// (d) MONOLITHIC 4-DOF coupled solve equals the one-way (STAGGERED) solve for a pure
// thermal-stress problem (task 4.1). No mechanical->thermal feedback (taylor_quinney
// == 0) => the 4-DOF tangent is block-triangular => the monolithic single-system solve
// reproduces the sequential result EXACTLY. This is the required monolithic validation.
void test_monolithic_equals_oneway() {
  Model m = coupled_bar_model();

  m.coupled_scheme = CoupledScheme::Staggered;
  const numerics::CoupledReport stag = numerics::solve_coupled_reported(m);
  m.coupled_scheme = CoupledScheme::Monolithic;
  const numerics::CoupledReport mono = numerics::solve_coupled_reported(m);

  // One-way: both converge in a single outer pass.
  CX_CHECK(stag.outer_iterations == 1);
  CX_CHECK(mono.outer_iterations == 1);
  // Temperature fields identical.
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i)
    CX_NEAR(mono.fields.thermal.temperature[i], stag.fields.thermal.temperature[i], 1e-9);
  // Displacement + stress fields identical to solver tolerance.
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i) {
    for (int c = 0; c < 3; ++c)
      CX_NEAR(mono.fields.mechanical.displacement[i][static_cast<std::size_t>(c)],
              stag.fields.mechanical.displacement[i][static_cast<std::size_t>(c)], 1e-9);
    for (int c = 0; c < 6; ++c)
      CX_NEAR(mono.fields.mechanical.stress[i][static_cast<std::size_t>(c)],
              stag.fields.mechanical.stress[i][static_cast<std::size_t>(c)], 1e-6);
  }
  // And both match the closed-form thermal stress on the hot face.
  const Real sxx_expect = -210000.0 * 1.2e-5 * 100.0;  // = -252.0
  Real rx_hi = 0.0;
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i)
    if (m.mesh.nodes()[i].x[0] == 1.0) rx_hi += mono.fields.mechanical.reaction[i][0];
  CX_NEAR(rx_hi, sxx_expect, 1e-2);
}

// (e) TWO-WAY plastic-dissipation heating (tasks 4.1 + 4.2). The mechanical->thermal
// coupling block: a single element is stretched past yield by a prescribed
// displacement, so plastic work accumulates; a Taylor-Quinney fraction of that work is
// deposited back into the thermal field as a volumetric heat source, raising the
// temperature above the purely-conducted profile. The STAGGERED scheme must iterate
// more than the single one-way pass (it discovers the dissipation source only after the
// first mechanical solve), and the MONOLITHIC scheme must converge to the SAME joint
// (T,u) state. This exercises the mechanical->thermal term the one-way scheme cannot.
void test_two_way_plastic_dissipation() {
  // One C3D8 hex. x=0 face held (u_x=0), x=1 face pulled to u_x=+0.01 (1% strain,
  // well past a 200 MPa yield), lateral rigid-body removed. One thermal BC pins the
  // field so conduction is well-posed; the plastic dissipation adds heat on top.
  const std::string deck = R"(
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
1, 1, 2, 3, 4, 5, 6, 7, 8
*MATERIAL, NAME=STEEL
*ELASTIC
210000., 0.3
*CONDUCTIVITY
50.
*EXPANSION
1.2e-5
*PLASTIC
200., 0.
*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL
*STEP
*COUPLED TEMPERATURE-DISPLACEMENT
*BOUNDARY
1, 11, 11, 100.
*BOUNDARY
1, 1, 1, 0.
4, 1, 1, 0.
5, 1, 1, 0.
8, 1, 1, 0.
2, 1, 1, 0.01
3, 1, 1, 0.01
6, 1, 1, 0.01
7, 1, 1, 0.01
1, 2
1, 3
2, 3
*END STEP
)";
  Model m = io::parse_inp(deck);
  CX_CHECK(m.has_plasticity());
  m.taylor_quinney = 0.9;  // deposit 90% of the plastic work as heat

  m.coupled_scheme = CoupledScheme::Staggered;
  const numerics::CoupledReport stag = numerics::solve_coupled_reported(m);
  m.coupled_scheme = CoupledScheme::Monolithic;
  const numerics::CoupledReport mono = numerics::solve_coupled_reported(m);

  // Two-way feedback: the staggered scheme iterates more than the single one-way pass.
  CX_CHECK(stag.outer_iterations > 1);
  // The element plastified: its mean equivalent plastic strain is > 0, so a positive
  // heat source was deposited and the temperature rose above the pinned 100 baseline.
  Real tmax = 0.0;
  for (const Real t : stag.fields.thermal.temperature) tmax = std::max(tmax, t);
  CX_CHECK(tmax > 100.0 + 1e-6);
  // Monolithic and staggered converge to the SAME joint temperature + displacement.
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i) {
    CX_NEAR(mono.fields.thermal.temperature[i], stag.fields.thermal.temperature[i], 1e-4);
    for (int c = 0; c < 3; ++c)
      CX_NEAR(mono.fields.mechanical.displacement[i][static_cast<std::size_t>(c)],
              stag.fields.mechanical.displacement[i][static_cast<std::size_t>(c)], 1e-7);
  }
}

}  // namespace

int main() {
  test_free_bar_expands_stress_free();
  test_constrained_bar_thermal_stress();
  test_constrained_bar_temp_dependent_alpha();
  test_no_temperature_is_pure_mechanics();
  test_coupled_deck_end_to_end();
  test_monolithic_equals_oneway();
  test_two_way_plastic_dissipation();
  if (cxtest::g_failures == 0) std::printf("test_thermal_stress: OK\n");
  CX_MAIN_RETURN();
}
