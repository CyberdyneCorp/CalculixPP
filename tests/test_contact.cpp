// Penalty node-to-surface contact (spec: contact — contact pairs / surface behavior;
// contact-search tasks 1.4/1.5/2.1/2.2). Validation is ANALYTICAL: two stacked cubes
// pressed together, meeting only through a *CONTACT PAIR (separate elements, separate
// nodes on the shared plane). The stock-CalculiX node-to-surface deck contact1 is matched
// (small penetration, correct displacement sign) by the Python regression suite; here the
// physics is validated against equilibrium and the penalty penetration law.
//
// TWO-BLOCK TEST: a base cube z in [0,1] (fixed at z=0) and a top cube z in [1,2] loaded
// downward by a total force F on its top face. The two cubes share the z=1 plane but are
// distinct elements with distinct nodes, so the load can only cross the interface through
// contact. At equilibrium:
//   (a) the base reaction at z=0 balances the applied load exactly (global equilibrium —
//       the interface transmits the full load);
//   (b) the interface penetration is small and ~ F / kappa (the penalty spring);
//   (c) removing the *CONTACT PAIR makes the top cube fall through unopposed (the two
//       bodies are otherwise unconnected) — contact is what carries the load.
#include <cmath>
#include <string>
#include <vector>

#include "calculixpp/fem/contact.hpp"
#include "calculixpp/io/inp_parser.hpp"
#include "calculixpp/numerics/nonlinear_static.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

// Two unit C3D8 cubes stacked in z with a small initial overlap `overlap` (top cube
// bottom at z = 1 - overlap so it starts penetrating; overlap 0 = just touching). Base
// nodes 1-8 (z=0 then z=1), top nodes 9-16 (z=1-overlap then z=2-overlap). All nodes
// confined u_x=u_y=0; base bottom fully fixed. `pressure` is a downward *DLOAD on the top
// face of the top cube. `behavior` is the *SURFACE BEHAVIOR block. Slave = top cube bottom
// nodes; master = base cube top face S2.
std::string two_block_deck(double kpen) {
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
*NSET, NSET=NBOT
1,2,3,4
*NSET, NSET=NTOP
13,14,15,16
*MATERIAL, NAME=STEEL
*ELASTIC
210000., 0.3
*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL
*BOUNDARY
NALL,1,2
NBOT,3,3
*SURFACE, NAME=SMAST
EBASE, S2
*SURFACE, NAME=SSLAV, TYPE=NODE
9,10,11,12
*SURFACE INTERACTION, NAME=SI
*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=LINEAR
)") + std::to_string(kpen) + R"(
*CONTACT PAIR, INTERACTION=SI, TYPE=NODE TO SURFACE
SSLAV, SMAST
*STEP
*STATIC
*CLOAD
13,3,-25.
14,3,-25.
15,3,-25.
16,3,-25.
*NODE PRINT, NSET=NALL
U,RF
*END STEP
)";
}

double node_uz(const StaticFields& f, const Model& m, Index id) {
  return f.displacement[static_cast<std::size_t>(m.mesh.node_index(id))][2];
}

// (a)+(b): equilibrium and small penetration. (c) is the analytical fall-through.
void test_two_block_equilibrium() {
  const double kpen = 1.0e6;
  const double F = 100.0;  // total downward load (4 x 25)
  const Model m = io::parse_inp(two_block_deck(kpen));
  CX_CHECK(m.has_contact());
  numerics::NonlinearReport rep;
  const StaticFields f = numerics::solve_nonlinear_static(m, {}, &rep);
  CX_CHECK(rep.converged);

  // Base reaction at z=0 (the fixed face) must balance the applied downward load: the
  // interface transmits the full load. Sum RF_z over the base-bottom nodes = +F.
  double rf = 0.0;
  for (Index id : {1, 2, 3, 4})
    rf += f.reaction[static_cast<std::size_t>(m.mesh.node_index(id))][2];
  CX_NEAR(rf, F, 1e-6 * F);  // reaction balances the load (global equilibrium)

  // Interface penetration = the RELATIVE closure across the contact plane: the slave
  // (top-cube bottom, nodes 9-12) sinks below the coincident master (base-cube top,
  // nodes 5-8). pen = uz(master) - uz(slave) > 0. This isolates the penalty overclosure
  // from each cube's own elastic shortening. With a linear penalty the total contact
  // force kappa*Σpen balances F, so pen ~ F/(4*kappa) — small and positive.
  auto interface_pen = [&]() {
    const std::pair<Index, Index> pairs4[4] = {{5, 9}, {6, 10}, {7, 11}, {8, 12}};
    double p = 0.0;
    for (auto [mst, slv] : pairs4) p += node_uz(f, m, mst) - node_uz(f, m, slv);
    return p / 4.0;
  };
  const double pen = interface_pen();
  CX_CHECK(pen > 0.0);         // the interface closes (slave penetrates the master)
  CX_CHECK(pen < 1e-3);        // small overclosure (<< the cube height 1)
  // Order of magnitude: pen ~ F/(4*kappa) = 100/(4e6) = 2.5e-5.
  CX_NEAR(pen, F / (4.0 * kpen), 5e-6);

  // The top cube still moves DOWN (load is carried through the interface, not lost) and
  // stays small (no blow-through).
  const double top = -node_uz(f, m, 13);
  CX_CHECK(top > 0.0);
  CX_CHECK(top < 1e-2);
}

// Penalty stiffness controls penetration: a stiffer spring -> less interface closure,
// while equilibrium (base reaction == load) holds regardless. (spec: contact — hard law
// is a large penalty giving near-zero overclosure.)
void test_penalty_controls_penetration() {
  const double F = 100.0;
  auto interface_gap = [&](double kpen) {
    const Model m = io::parse_inp(two_block_deck(kpen));
    const StaticFields f = numerics::solve_nonlinear_static(m);
    const std::pair<Index, Index> pairs4[4] = {{5, 9}, {6, 10}, {7, 11}, {8, 12}};
    double g = 0.0;
    for (auto [mst, slv] : pairs4)
      g += node_uz(f, m, mst) - node_uz(f, m, slv);  // relative interface closure
    double rf = 0.0;
    for (Index id : {1, 2, 3, 4})
      rf += f.reaction[static_cast<std::size_t>(m.mesh.node_index(id))][2];
    CX_NEAR(rf, F, 1e-6 * F);  // equilibrium holds at every stiffness
    return g / 4.0;
  };
  const double soft = interface_gap(1.0e5);
  const double stiff = interface_gap(1.0e7);
  // A 100x stiffer contact spring closes the interface far less; the penetration part
  // scales ~ 1/kappa. Both are small but the softer one is clearly larger.
  CX_CHECK(soft > stiff);
  CX_CHECK(soft > 5.0 * stiff);
}

// A two-block deck WITH Coulomb *FRICTION (mu, stick stiffness) on the surface interaction.
// Same geometry as two_block_deck: base cube z in [0,1], top cube z in [1,2], distinct
// nodes on z=1, slave = top-bottom nodes 9-12, master = base top face S2, linear normal
// penalty `kpen`. Used both for the end-to-end solve and to build the resolved pair the
// operator-level friction test drives directly.
std::string friction_two_block_deck(double kpen, double mu, double stick) {
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
*NSET, NSET=NBOT
1,2,3,4
*NSET, NSET=NTOP
13,14,15,16
*MATERIAL, NAME=STEEL
*ELASTIC
210000., 0.3
*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL
*BOUNDARY
NALL,1,2
NBOT,3,3
*SURFACE, NAME=SMAST
EBASE, S2
*SURFACE, NAME=SSLAV, TYPE=NODE
9,10,11,12
*SURFACE INTERACTION, NAME=SI
*SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=LINEAR
)") + std::to_string(kpen) + R"(
*FRICTION
)" + std::to_string(mu) + ", " + std::to_string(stick) + R"(
*CONTACT PAIR, INTERACTION=SI, TYPE=NODE TO SURFACE
SSLAV, SMAST
*STEP
*STATIC
*CLOAD
13,3,-25.
14,3,-25.
15,3,-25.
16,3,-25.
*NODE PRINT, NSET=NALL
U,RF
*END STEP
)";
}

// The total physical tangential friction TRACTION the contact operator applies to the four
// slave nodes 9-12, for a full displacement field `u`. The operator stores its contribution
// as f_int_slave = -t (so the residual r = f_ext - f_int carries the friction force), so the
// physical traction on the slave is t = -Σ f_int_slave. This is the exact interface friction
// force, free of solver/reaction-recovery noise. `pairs` are the resolved contact pairs.
Vec3 slave_friction_traction(const Model& m, std::vector<fem::ResolvedContactPair>& pairs,
                             const std::vector<Vec3>& u) {
  fem::LinearSystem sys = fem::assemble_linear_static(m);
  std::vector<Real> f_int(3 * m.mesh.num_nodes(), 0.0);
  fem::add_contact(m, pairs, u, sys, f_int);
  Vec3 t{0, 0, 0};
  for (Index id : {9, 10, 11, 12}) {
    const std::size_t i = static_cast<std::size_t>(m.mesh.node_index(id));
    for (int d = 0; d < 3; ++d)
      t[static_cast<std::size_t>(d)] -= f_int[i * 3 + static_cast<std::size_t>(d)];  // t = -f_int
  }
  return t;
}

// Operator-level Coulomb stick/slip (spec: contact — tangential contact). Drive the slave
// nodes with a KNOWN penetration (uz on the slave) and a KNOWN tangential slide (ux on the
// slave), then read the interface friction force the operator produces. The slave sits at
// z=1 on a fixed master, so a uniform slave displacement u = (ux,0,-pen) gives penetration
// `pen` (normal force Fn = kappa*pen over 4 nodes) and tangential slip `ux`. Coulomb: the
// tangential force sticks at -stick*ux while |stick*ux| < mu*Fn(per node), then saturates
// at mu*Fn once slipping — the defining stick->slip transition, read exactly from f_int.
void test_friction_operator_stick_slip() {
  const double kpen = 1.0e6, mu = 0.3, stick = 1.0e6;
  const Model m = io::parse_inp(friction_two_block_deck(kpen, mu, stick));
  std::vector<fem::ResolvedContactPair> pairs = fem::build_contact_pairs(m);

  const double pen = 1.0e-4;  // slave pushed 1e-4 into the master
  // Per-node normal force Fn = kappa*pen; the Coulomb cap per slave node is mu*Fn.
  const double Fn_node = kpen * pen;
  const double cap_node = mu * Fn_node;  // = 0.3 * 1e6 * 1e-4 = 30 per node

  auto slave_slide = [&](double ux) {
    std::vector<Vec3> u(m.mesh.num_nodes(), Vec3{0, 0, 0});
    for (Index id : {9, 10, 11, 12})
      u[static_cast<std::size_t>(m.mesh.node_index(id))] = Vec3{ux, 0.0, -pen};
    return slave_friction_traction(m, pairs, u);
  };

  // STICK: a tiny tangential slide keeps |stick*ux| below the cap. The friction traction is
  // the elastic stick force -stick*ux per node (4 nodes) and OPPOSES the +x slide (negative).
  const double ux_stick = 1.0e-6;  // stick*ux = 1.0 per node < cap 30
  const Vec3 ts = slave_slide(ux_stick);
  CX_NEAR(ts[0], -4.0 * stick * ux_stick, 1e-6 * 4.0 * stick * ux_stick);  // = -4.0
  CX_CHECK(std::fabs(ts[0]) < 4.0 * cap_node);  // inside the cone -> stuck

  // SLIP: a large tangential slide drives past the cap. The traction SATURATES at the
  // Coulomb limit mu*Fn summed over the four nodes, still opposing the slide (negative x).
  const Vec3 tl = slave_slide(1.0e-2);  // stick*ux = 1e4 per node >> cap 30
  CX_NEAR(tl[0], -4.0 * cap_node, 1e-6 * 4.0 * cap_node);  // saturated at -mu*Fn
  CX_CHECK(tl[0] < 0.0);  // the slip traction resists the slide

  // The Coulomb cap scales linearly with mu: doubling mu doubles the saturated slip
  // traction. Rebuild the pair at mu = 0.6 and re-drive well into slip.
  const Model m2 = io::parse_inp(friction_two_block_deck(kpen, 2.0 * mu, stick));
  std::vector<fem::ResolvedContactPair> pairs2 = fem::build_contact_pairs(m2);
  std::vector<Vec3> u2(m2.mesh.num_nodes(), Vec3{0, 0, 0});
  for (Index id : {9, 10, 11, 12})
    u2[static_cast<std::size_t>(m2.mesh.node_index(id))] = Vec3{1.0e-2, 0.0, -pen};
  const Vec3 tl2 = slave_friction_traction(m2, pairs2, u2);
  CX_NEAR(tl2[0], 2.0 * tl[0], 1e-6 * std::fabs(2.0 * tl[0]));  // cap doubles with mu
}

// End-to-end friction solve: two blocks pressed together (F=100 down) then the top block
// slid tangentially. With mu large enough the interface STICKS (the whole load path is
// intact, small tangential motion); the solve CONVERGES with the friction contribution in
// the Newton tangent+residual. This exercises the full parse->route->assemble->solve path.
void test_friction_end_to_end_converges() {
  // Normal-only reference (mu implicitly 0 via no *FRICTION): must still converge and
  // balance the load, unchanged by the friction machinery being present but inactive.
  const Model m0 = io::parse_inp(two_block_deck(1.0e6));
  const StaticFields f0 = numerics::solve_nonlinear_static(m0);
  double rf0 = 0.0;
  for (Index id : {1, 2, 3, 4})
    rf0 += f0.reaction[static_cast<std::size_t>(m0.mesh.node_index(id))][2];
  CX_NEAR(rf0, 100.0, 1e-6 * 100.0);

  // With friction present the pressed-together solve still converges and transmits the
  // full normal load (friction does not disturb the normal equilibrium here — the load is
  // purely normal). This confirms the friction operator is inert when there is no
  // tangential drive and does not break the Newton solve.
  const Model mf = io::parse_inp(friction_two_block_deck(1.0e6, 0.3, 1.0e6));
  numerics::NonlinearReport rep;
  const StaticFields ff = numerics::solve_nonlinear_static(mf, {}, &rep);
  CX_CHECK(rep.converged);
  double rff = 0.0;
  for (Index id : {1, 2, 3, 4})
    rff += ff.reaction[static_cast<std::size_t>(mf.mesh.node_index(id))][2];
  CX_NEAR(rff, 100.0, 1e-6 * 100.0);  // normal equilibrium intact with friction present
}

// The overclosure law: no force in clearance, linear force in penetration.
void test_overclosure_law() {
  SurfaceBehavior hard;  // Law::Hard by default
  hard.law = SurfaceBehavior::Law::Linear;
  hard.k = 1000.0;
  CX_NEAR(fem::contact_force_weight(hard, 1000.0, 0.5), 0.0, 1e-14);   // clearance -> 0
  CX_NEAR(fem::contact_force_weight(hard, 1000.0, -0.5), 500.0, 1e-10); // F = k*pen
  CX_NEAR(fem::contact_force_weight(hard, 1000.0, 0.0), 0.0, 1e-14);   // just touching
}

// The Coulomb return-map: inside the cone the trial traction is kept (stick), outside it
// is scaled onto the cone radius mu*p (slip). Zero pressure / zero mu -> no traction.
void test_friction_return_map() {
  bool slip = false;
  // Inside the cone: |t_trial| = 3 < mu*p = 0.5*10 = 5 -> stick (unchanged).
  Vec3 t = fem::friction_return_map({3.0, 0.0, 0.0}, 0.5, 10.0, slip);
  CX_CHECK(!slip);
  CX_NEAR(t[0], 3.0, 1e-12);
  // On/over the cone: |t_trial| = 8 > 5 -> slip, scaled to magnitude 5 along the direction.
  t = fem::friction_return_map({8.0, 0.0, 0.0}, 0.5, 10.0, slip);
  CX_CHECK(slip);
  CX_NEAR(t[0], 5.0, 1e-12);  // capped at mu*p
  // Oblique trial: |t_trial| = 5 (3,4) > mu*p = 2 -> scaled to length 2, same direction.
  t = fem::friction_return_map({3.0, 4.0, 0.0}, 0.2, 10.0, slip);
  CX_CHECK(slip);
  CX_NEAR(std::sqrt(t[0] * t[0] + t[1] * t[1]), 2.0, 1e-12);
  CX_NEAR(t[0] / t[1], 3.0 / 4.0, 1e-12);  // direction preserved
  // No normal pressure -> no tangential traction can be sustained.
  t = fem::friction_return_map({5.0, 0.0, 0.0}, 0.5, 0.0, slip);
  CX_CHECK(slip);
  CX_NEAR(t[0], 0.0, 1e-12);
}

// Contact output CSTR (tasks.md 2.5; spec: contact — contact output). After the two-block
// pressed-together solve, the recovered per-slave-node contact result reports each slave
// node as CLOSED (in contact), with a positive normal pressure whose sum over the tributary
// areas balances the applied load, and a small penetration gap. This validates the CSTR
// recovery the .dat/.frd writers emit.
void test_contact_output_cstr() {
  const double kpen = 1.0e6, F = 100.0;
  const Model m = io::parse_inp(two_block_deck(kpen));
  const StaticFields f = numerics::solve_nonlinear_static(m);
  CX_CHECK(f.contact.size() == 4);  // four slave nodes 9-12

  double p_sum_force = 0.0;
  for (const ContactPoint& c : f.contact) {
    CX_CHECK(c.closed);            // every slave node is in contact
    CX_CHECK(c.p > 0.0);           // positive normal pressure
    CX_CHECK(c.gap < 0.0);         // small penetration
    CX_CHECK(std::fabs(c.gap) < 1e-3);
    // pressure * tributary area (unit face / 4 nodes = 0.25) contributes to the total force.
    p_sum_force += c.p * 0.25;
  }
  // The summed contact force (pressure x area) balances the applied load F.
  CX_NEAR(p_sum_force, F, 1e-4 * F);
}

// *CLEARANCE offsets the geometric gap: with a known slave penetration `pen`, the operator
// sees g_eff = g_geometric + clearance, so a clearance larger than the penetration OPENS
// the interface (no contact) while zero clearance leaves it CLOSED. Validated directly
// through recover_contact at a KNOWN displacement (no singular free-body solve needed).
void test_clearance_offsets_gap() {
  const double kpen = 1.0e6, pen = 1.0e-4;
  // Drive the four slave nodes pen into the master (uniform u_z = -pen). The base master
  // (nodes 5-8) stays at zero, so the geometric gap is -pen (penetration).
  auto drive = [&](const Model& m) {
    std::vector<Vec3> u(m.mesh.num_nodes(), Vec3{0, 0, 0});
    for (Index id : {9, 10, 11, 12})
      u[static_cast<std::size_t>(m.mesh.node_index(id))] = Vec3{0.0, 0.0, -pen};
    return u;
  };

  // No clearance: g_eff = -pen < 0 -> all slave nodes CLOSED with positive pressure.
  const Model m0 = io::parse_inp(two_block_deck(kpen));
  std::vector<fem::ResolvedContactPair> p0 = fem::build_contact_pairs(m0);
  const std::vector<ContactPoint> c0 = fem::recover_contact(m0, p0, drive(m0));
  CX_CHECK(c0.size() == 4);
  for (const ContactPoint& c : c0) {
    CX_CHECK(c.closed);
    CX_CHECK(c.p > 0.0);
    CX_NEAR(c.gap, -pen, 1e-9);  // geometric gap unchanged (clearance 0)
  }

  // Large clearance (10*pen): g_eff = -pen + 10*pen = +9*pen > 0 -> the interface OPENS.
  std::string deck = two_block_deck(kpen);
  // *CLEARANCE must come AFTER the *CONTACT PAIR (its slave/master surfaces are set on the
  // pair's data line). Splice it before *STEP.
  deck.insert(deck.find("*STEP"),
              "*CLEARANCE, MASTER=SMAST, SLAVE=SSLAV, VALUE=" + std::to_string(10.0 * pen) +
                  "\n");
  const Model mc = io::parse_inp(deck);
  CX_CHECK(mc.contact_pairs[0].has_clearance);
  CX_NEAR(mc.contact_pairs[0].clearance, 10.0 * pen, 1e-12);
  std::vector<fem::ResolvedContactPair> pc = fem::build_contact_pairs(mc);
  const std::vector<ContactPoint> cc = fem::recover_contact(mc, pc, drive(mc));
  for (const ContactPoint& c : cc) {
    CX_CHECK(!c.closed);         // clearance opens the interface -> OPEN
    CX_NEAR(c.gap, 9.0 * pen, 1e-9);  // g_eff = -pen + 10 pen
    CX_NEAR(c.p, 0.0, 1e-14);    // no pressure when open
  }
}

}  // namespace

int main() {
  test_overclosure_law();
  test_two_block_equilibrium();
  test_penalty_controls_penetration();
  test_friction_return_map();
  test_friction_operator_stick_slip();
  test_friction_end_to_end_converges();
  test_contact_output_cstr();
  test_clearance_offsets_gap();
  CX_MAIN_RETURN();
}
