// J2 (von Mises) plasticity tests (spec: material-models 4.1/4.2).
//   1. Below yield the model reduces EXACTLY to linear elasticity (stress + tangent).
//   2. Radial return lands the trial stress on the yield surface (von Mises == sy).
//   3. The returned 6x6 consistent tangent matches a finite-difference tangent
//      (perturb each strain component) — the 4.2 verification requirement.
//   4. Kinematic hardening moves the back stress (COMBINED / KINEMATIC).
//   5. Committed vs. trial state: evaluate() never advances committed history; the
//      driver's commit() promotes trial -> committed.
//   6. End-to-end single-hex uniaxial tension follows the analytic hardening curve
//      (validation (a): stress = sy0 + H_iso * eps_p, exactly known).
#include <array>
#include <cmath>
#include <vector>

#include "calculixpp/fem/assembly.hpp"
#include "calculixpp/fem/material_model.hpp"
#include "calculixpp/numerics/nonlinear_static.hpp"
#include "check.hpp"

using namespace cxpp;
using namespace cxpp::fem;

namespace {

ElasticIso steel() { return ElasticIso{210000.0, 0.3}; }

// Isotropic hardening: sy0=800 at ep=0, sy=960 at ep=0.02 -> H_iso = 8000.
Plastic iso_hardening() {
  Plastic p;
  p.hardening = Plastic::Hardening::Isotropic;
  p.eqplastic = {0.0, 0.02};
  p.yield = {800.0, 960.0};
  return p;
}

Real von_mises(const Voigt6& s) {
  const Real p = (s[0] + s[1] + s[2]) / 3.0;
  const Real dxx = s[0] - p, dyy = s[1] - p, dzz = s[2] - p;
  const Real j2 = 0.5 * (dxx * dxx + dyy * dyy + dzz * dzz) +
                  s[3] * s[3] + s[4] * s[4] + s[5] * s[5];
  return std::sqrt(3.0 * j2);
}

void test_elastic_below_yield() {
  // A tiny strain stays inside the yield surface: stress and tangent must equal the
  // linear-elastic model bit-for-bit (elastic decks are unaffected by plasticity).
  J2PlasticMaterial jp(steel(), iso_hardening());
  ElasticIsoMaterial el(steel());
  const Voigt6 strain{1e-4, -3e-5, -3e-5, 2e-5, 0.0, 1e-5};
  MaterialState s1, s2;
  const MaterialResponse rp = jp.evaluate(strain, s1);
  const MaterialResponse re = el.evaluate(strain, s2);
  for (int i = 0; i < 6; ++i) {
    CX_NEAR(rp.stress[static_cast<std::size_t>(i)],
            re.stress[static_cast<std::size_t>(i)], 1e-9);
    for (int j = 0; j < 6; ++j)
      CX_NEAR(rp.tangent[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)],
              re.tangent[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)],
              1e-9);
  }
  CX_CHECK(s1.eqplastic == 0.0);  // no plastic flow
}

void test_return_to_yield_surface() {
  // Perfect plasticity: a strain well past yield must return the stress exactly onto
  // the von Mises surface at the (constant) yield stress.
  Plastic perf;
  perf.eqplastic = {0.0};
  perf.yield = {800.0};
  J2PlasticMaterial jp(steel(), perf);
  const Voigt6 strain{0.01, -0.003, -0.003, 0.0, 0.0, 0.0};
  MaterialState st;
  const MaterialResponse r = jp.evaluate(strain, st);
  CX_NEAR(von_mises(r.stress), 800.0, 1e-6);
  CX_CHECK(st.eqplastic > 0.0);

  // Isotropic hardening: von Mises must equal the grown yield sy0 + H_iso * eqp.
  J2PlasticMaterial jh(steel(), iso_hardening());
  MaterialState st2;
  const MaterialResponse rh = jh.evaluate(strain, st2);
  const Real sy_expected = 800.0 + 8000.0 * st2.eqplastic;
  CX_NEAR(von_mises(rh.stress), sy_expected, 1e-5);
}

// Finite-difference tangent verification (4.2): perturb each strain component about a
// plastic state (from the SAME committed history each time) and compare the numerical
// d sigma / d strain to the returned consistent tangent.
void test_consistent_tangent_fd() {
  const std::array<Plastic, 2> models{
      iso_hardening(), [] {
        Plastic p = iso_hardening();
        p.hardening = Plastic::Hardening::Combined;
        p.kinematic_modulus = 5000.0;
        return p;
      }()};
  for (const Plastic& pl : models) {
    J2PlasticMaterial jp(steel(), pl);
    const Voigt6 strain{0.008, -0.0015, -0.002, 0.0012, 0.0006, -0.0004};
    MaterialState base;  // virgin committed history
    const MaterialResponse r = jp.evaluate(strain, base);
    CX_CHECK(base.eqplastic > 0.0);  // ensure we are on the plastic branch

    const Real h = 1e-7;
    for (int j = 0; j < 6; ++j) {
      Voigt6 ep = strain, em = strain;
      ep[static_cast<std::size_t>(j)] += h;
      em[static_cast<std::size_t>(j)] -= h;
      MaterialState cp, cm;  // both virgin -> same committed state
      const Voigt6 sp = jp.evaluate(ep, cp).stress;
      const Voigt6 sm = jp.evaluate(em, cm).stress;
      for (int i = 0; i < 6; ++i) {
        const Real fd = (sp[static_cast<std::size_t>(i)] -
                         sm[static_cast<std::size_t>(i)]) / (2 * h);
        const Real an = r.tangent[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
        CX_CHECK(std::abs(fd - an) <= 1e-3 * (std::abs(an) + 1.0));
      }
    }
  }
}

void test_kinematic_back_stress() {
  // Pure kinematic hardening: the yield surface does not grow (isotropic yield stays
  // at sy0), but the back stress moves, so the relative stress magnitude stays at sy0.
  Plastic pk;
  pk.hardening = Plastic::Hardening::Kinematic;
  pk.eqplastic = {0.0, 0.02};
  pk.yield = {800.0, 960.0};  // slope 8000 -> H_kin
  J2PlasticMaterial jp(steel(), pk);
  const Voigt6 strain{0.01, -0.003, -0.003, 0.0, 0.0, 0.0};
  MaterialState st;
  const MaterialResponse r = jp.evaluate(strain, st);
  // Back stress advanced; von Mises of (dev(sigma) - alpha) equals sy0.
  bool alpha_nonzero = false;
  for (int i = 0; i < 6; ++i)
    if (std::abs(st.back_stress[static_cast<std::size_t>(i)]) > 1e-6) alpha_nonzero = true;
  CX_CHECK(alpha_nonzero);
  // relative stress xi = dev(sigma) - alpha; its von Mises must be sy0 (surface fixed)
  const Real p = (r.stress[0] + r.stress[1] + r.stress[2]) / 3.0;
  Voigt6 xi{r.stress[0] - p - st.back_stress[0], r.stress[1] - p - st.back_stress[1],
            r.stress[2] - p - st.back_stress[2], r.stress[3] - st.back_stress[3],
            r.stress[4] - st.back_stress[4], r.stress[5] - st.back_stress[5]};
  CX_NEAR(von_mises(xi), 800.0, 1e-4);
}

void test_commit_semantics() {
  // evaluate() writes trial history, never committed; commit() promotes it.
  J2PlasticMaterial jp(steel(), iso_hardening());
  const Voigt6 strain{0.01, -0.003, -0.003, 0.0, 0.0, 0.0};
  MaterialState st;
  jp.evaluate(strain, st);
  CX_CHECK(st.eqplastic > 0.0);
  CX_CHECK(st.committed_eqplastic == 0.0);  // committed untouched by evaluate
  st.commit();
  CX_NEAR(st.committed_eqplastic, st.eqplastic, 1e-15);
  // Re-evaluating the same strain from the now-committed state yields ~zero extra
  // plastic flow (already on the surface for that strain).
  const Real ep_before = st.eqplastic;
  jp.evaluate(strain, st);
  CX_CHECK(st.eqplastic >= ep_before - 1e-12);
}

// End-to-end validation (a): a single C3D8 cube in uniaxial tension. Symmetry BCs on
// the three back faces, a prescribed axial displacement on the top; the recovered
// axial stress must follow the analytic uniaxial hardening curve
// sigma = sy0 + H_iso * eps_p with eps_p = (E*eps - sy0)/(E + H_iso).
Model uniaxial_cube(Real uz) {
  Model m;
  const std::array<Vec3, 8> x{{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                               {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}};
  for (int i = 0; i < 8; ++i) m.mesh.add_node(i + 1, x[static_cast<std::size_t>(i)]);
  m.mesh.add_element(1, ElementType::C3D8, {1, 2, 3, 4, 5, 6, 7, 8});
  m.mesh.add_elset("EALL", {1});
  Material mat{"STEEL", steel(), std::nullopt, iso_hardening()};
  m.materials["STEEL"] = mat;
  m.sections.push_back(SolidSection{"EALL", "STEEL"});
  for (Index n : {1, 4, 5, 8}) m.spcs.push_back(Spc{n, 1, 0.0});  // x=0 face ux=0
  for (Index n : {1, 2, 5, 6}) m.spcs.push_back(Spc{n, 2, 0.0});  // y=0 face uy=0
  for (Index n : {1, 2, 3, 4}) m.spcs.push_back(Spc{n, 3, 0.0});  // z=0 face uz=0
  for (Index n : {5, 6, 7, 8}) m.spcs.push_back(Spc{n, 3, uz});   // top face uz
  return m;
}

void test_uniaxial_end_to_end() {
  const Real E = 210000.0, H = 8000.0, sy0 = 800.0;
  const Real uz = 0.01;  // axial strain 0.01, well past yield (eps_y = sy0/E ~ 0.0038)
  const Model m = uniaxial_cube(uz);
  numerics::NonlinearReport rep;
  const StaticFields f = numerics::solve_nonlinear_static(m, {}, &rep);
  CX_CHECK(rep.converged);

  Real szz = 0.0;
  for (const Voigt6& s : f.stress) szz += s[2];
  szz /= static_cast<Real>(f.stress.size());

  const Real eps_p = (E * uz - sy0) / (E + H);
  const Real sigma = sy0 + H * eps_p;
  CX_CHECK(eps_p > 0.0);
  CX_NEAR(szz, sigma, 1e-4 * sigma);
}

}  // namespace

int main() {
  test_elastic_below_yield();
  test_return_to_yield_surface();
  test_consistent_tangent_fd();
  test_kinematic_back_stress();
  test_commit_semantics();
  test_uniaxial_end_to_end();
  if (cxtest::g_failures == 0) std::printf("test_plasticity: OK\n");
  CX_MAIN_RETURN();
}
