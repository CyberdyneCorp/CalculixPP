// Stress-based high-cycle-fatigue (HCF) driver (add-high-cycle-fatigue). Covers:
//   * the analytical stress-life oracle: for a point of recovered stress amplitude S and a
//     Basquin S-N curve (a, b), the reported life equals the closed form N = (S/a)^(1/b) to
//     < 1e-6 relative, and the worst-case node is the highest-amplitude node;
//   * the criterion selection (signed vs plain von Mises);
//   * the missing-source diagnostics (no *FATIGUE curve; no elastic material).
#include <cmath>
#include <stdexcept>
#include <vector>

#include "calculixpp/io/inp_parser.hpp"
#include "calculixpp/numerics/high_cycle_fatigue.hpp"
#include "calculixpp/numerics/linear_static.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

const std::vector<Vec3> kUnitCube = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                     {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

// von Mises equivalent (mirrors the driver's reduction) for the oracle.
Real von_mises(const Voigt6& s) {
  const Real a = s[0] - s[1], b = s[1] - s[2], c = s[2] - s[0];
  return std::sqrt(0.5 * (a * a + b * b + c * c) +
                   3.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]));
}

// A uniaxial-tension C3D8 cube: fixed at x=0 face (roller in x on nodes 1,4,5,8), pulled by
// a uniform +x traction on the x=1 face (nodes 2,3,6,7 loaded in +x). One node pinned in
// y/z to remove rigid modes. This yields a (near-)uniform uniaxial σxx tension state — the
// stress amplitude the HCF driver reduces from.
Model make_uniaxial_cube(const SNCurve* sn, bool with_elastic = true) {
  Model m;
  for (std::size_t i = 0; i < kUnitCube.size(); ++i)
    m.mesh.add_node(static_cast<Index>(i + 1), kUnitCube[i]);
  m.mesh.add_element(1, ElementType::C3D8, {1, 2, 3, 4, 5, 6, 7, 8});
  m.mesh.add_elset("EALL", {1});

  Material mat;
  mat.name = "STEEL";
  if (with_elastic) mat.elastic = ElasticIso{210000.0, 0.3};
  if (sn) mat.sn_curve = *sn;
  m.materials["STEEL"] = mat;
  m.sections.push_back(SolidSection{"EALL", "STEEL"});

  // x=0 face fixed in x; pin node 1 fully to kill rigid-body y/z translation/rotation.
  for (Index nd : {1, 4, 5, 8}) m.spcs.push_back(Spc{nd, 1, 0.0});
  m.spcs.push_back(Spc{1, 2, 0.0});
  m.spcs.push_back(Spc{1, 3, 0.0});
  m.spcs.push_back(Spc{4, 3, 0.0});
  m.spcs.push_back(Spc{5, 2, 0.0});
  // +x traction on the x=1 face: total force F spread over 4 nodes (unit face area).
  const Real F = 100.0;  // total axial force -> nominal σxx = F / A = 100 (A = 1)
  for (Index nd : {2, 3, 6, 7}) m.cloads.push_back(Cload{nd, 1, F / 4.0});

  m.procedure = Procedure::HighCycleFatigue;
  return m;
}

// Analytical oracle: the reported worst-case life equals the closed-form Basquin inversion
// N = (S/a)^(1/b) of the recovered worst-case amplitude, and the worst node is the
// highest-amplitude node.
void test_analytical_life() {
  const SNCurve sn{1000.0, -0.1};  // S_a = 1000 * N^(-0.1)
  const Model m = make_uniaxial_cube(&sn);
  const numerics::HcfReport rep = numerics::evaluate_hcf(m);

  // Independently recover the stress field and find the max-amplitude node.
  const StaticFields f = numerics::solve_linear_static(m);
  Real max_amp = 0.0;
  std::size_t max_i = 0;
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i) {
    // Signed von Mises (default criterion): tension -> positive mean -> +vm.
    const Real vm = von_mises(f.stress[i]);
    if (vm > max_amp) {
      max_amp = vm;
      max_i = i;
    }
  }
  CX_CHECK(max_amp > 0.0);

  // Worst node is the highest-amplitude node.
  CX_CHECK(rep.worst_index == max_i);
  CX_CHECK(rep.worst_node_id == m.mesh.nodes()[max_i].id);
  CX_NEAR(rep.worst_amplitude, max_amp, 1e-9 * max_amp);

  // Closed-form Basquin life at the recovered amplitude, matched to < 1e-6 relative.
  const Real n_closed = std::pow(max_amp / sn.a, 1.0 / sn.b);
  CX_NEAR(rep.worst_life, n_closed, 1e-6 * n_closed);

  // And every node's life matches its own closed-form inversion.
  for (std::size_t i = 0; i < m.mesh.num_nodes(); ++i) {
    const Real s_a = std::fabs(rep.amplitude[i]);
    if (s_a <= 0.0) continue;
    const Real n_i = std::pow(s_a / sn.a, 1.0 / sn.b);
    CX_NEAR(rep.life[i], n_i, 1e-6 * n_i);
  }
}

// A larger amplitude gives a shorter life (b < 0 monotonicity): scaling the load up must
// reduce the reported worst-case life.
void test_amplitude_monotonicity() {
  const SNCurve sn{1000.0, -0.1};
  const numerics::HcfReport lo = numerics::evaluate_hcf(make_uniaxial_cube(&sn));

  Model hi = make_uniaxial_cube(&sn);
  for (Cload& cl : hi.cloads) cl.value *= 2.0;  // double the axial load
  const numerics::HcfReport hr = numerics::evaluate_hcf(hi);

  CX_CHECK(hr.worst_amplitude > lo.worst_amplitude);
  CX_CHECK(hr.worst_life < lo.worst_life);
}

// Missing S-N curve -> diagnostic error, no output.
void test_missing_sn_curve() {
  const Model m = make_uniaxial_cube(nullptr);  // no *FATIGUE curve
  bool threw = false;
  try {
    numerics::evaluate_hcf(m);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CX_CHECK(threw);
}

// Missing elastic material (no stress source) -> diagnostic error.
void test_missing_stress_source() {
  const SNCurve sn{1000.0, -0.1};
  const Model m = make_uniaxial_cube(&sn, /*with_elastic=*/false);
  bool threw = false;
  try {
    numerics::evaluate_hcf(m);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  CX_CHECK(threw);
}

// End-to-end through the parser: a *HCF + *FATIGUE deck parses to the HCF procedure and
// evaluates to a finite worst-case life.
void test_parse_and_evaluate() {
  const char* deck = R"(
*NODE
1, 0.0, 0.0, 0.0
2, 1.0, 0.0, 0.0
3, 1.0, 1.0, 0.0
4, 0.0, 1.0, 0.0
5, 0.0, 0.0, 1.0
6, 1.0, 0.0, 1.0
7, 1.0, 1.0, 1.0
8, 0.0, 1.0, 1.0
*ELEMENT, TYPE=C3D8, ELSET=EALL
1, 1, 2, 3, 4, 5, 6, 7, 8
*MATERIAL, NAME=STEEL
*ELASTIC
210000.0, 0.3
*FATIGUE
1000.0, -0.1
*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL
*STEP
*HCF, CRITERION=SIGNED-VON-MISES
*BOUNDARY
1, 1, 3
4, 1, 1
4, 3, 3
5, 1, 1
5, 2, 2
8, 1, 1
*CLOAD
2, 1, 25.0
3, 1, 25.0
6, 1, 25.0
7, 1, 25.0
*END STEP
)";
  const Model m = io::parse_inp(deck);
  CX_CHECK(m.procedure == Procedure::HighCycleFatigue);
  CX_CHECK(m.hcf_criterion == FatigueCriterion::SignedVonMises);
  const numerics::HcfReport rep = numerics::evaluate_hcf(m);
  CX_CHECK(rep.worst_life > 0.0);
  CX_CHECK(std::isfinite(rep.worst_life));
  CX_CHECK(rep.worst_amplitude > 0.0);
}

}  // namespace

int main() {
  test_analytical_life();
  test_amplitude_monotonicity();
  test_missing_sn_curve();
  test_missing_stress_source();
  test_parse_and_evaluate();
  CX_MAIN_RETURN();
}
