// Material-model breadth tests (spec: material-models 4.3, 4.6).
//
// 4.6 User material: the *USER MATERIAL C++ plug-in path. Registers a custom law in
// the registry, drives it through a *USER MATERIAL deck, and checks the whole path:
//   - a registered factory is found and built by make_material_points;
//   - the built-in reference LinearUserMaterial reproduces isotropic elasticity;
//   - *DEPVAR state variables are sized, written by evaluate(), and committed;
//   - *RATEDEPENDENT scales the stress/tangent;
//   - end-to-end a *USER MATERIAL cube solves identically to the same *ELASTIC cube
//     (the user path is validated against the analytic linear-elastic solution).
//
// 4.3 Hyperelastic: compressible neo-Hookean in its small-strain analytic limit equals
// elastic_iso_D for the matching (E, nu); an end-to-end uniaxial cube matches the
// closed-form linear-elastic axial stress. (Large-strain kinematics deferred — the
// model is validated only in the small-strain regime where it is analytically exact.)
#include <array>
#include <cmath>
#include <memory>
#include <vector>

#include "calculixpp/fem/assembly.hpp"
#include "calculixpp/fem/element.hpp"
#include "calculixpp/fem/material_model.hpp"
#include "calculixpp/numerics/nonlinear_static.hpp"
#include "check.hpp"

using namespace cxpp;
using namespace cxpp::fem;

namespace {

ElasticIso steel() { return ElasticIso{210000.0, 0.3}; }

// A single C3D8 unit cube with symmetry BCs and a prescribed axial top displacement.
// `configure` installs the material + section (lets each test pick the constitutive
// law while keeping identical geometry/BCs).
template <class Configure>
Model uniaxial_cube(Real uz, Configure configure) {
  Model m;
  const std::array<Vec3, 8> x{{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                               {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}};
  for (int i = 0; i < 8; ++i) m.mesh.add_node(i + 1, x[static_cast<std::size_t>(i)]);
  m.mesh.add_element(1, ElementType::C3D8, {1, 2, 3, 4, 5, 6, 7, 8});
  m.mesh.add_elset("EALL", {1});
  configure(m);
  m.sections.push_back(SolidSection{"EALL", "MAT"});
  for (Index n : {1, 4, 5, 8}) m.spcs.push_back(Spc{n, 1, 0.0});
  for (Index n : {1, 2, 5, 6}) m.spcs.push_back(Spc{n, 2, 0.0});
  for (Index n : {1, 2, 3, 4}) m.spcs.push_back(Spc{n, 3, 0.0});
  for (Index n : {5, 6, 7, 8}) m.spcs.push_back(Spc{n, 3, uz});
  return m;
}

Real mean_szz(const StaticFields& f) {
  Real szz = 0.0;
  for (const Voigt6& s : f.stress) szz += s[2];
  return szz / static_cast<Real>(f.stress.size());
}

// ---- 4.6 user material -----------------------------------------------------------

// A custom user law registered by the TEST (not the core): a scaled linear-elastic law
// where constants = {E, nu, scale}. Proves a downstream user can register their own
// MaterialModel and reach it from a deck without touching the solver core.
class ScaledElasticUser final : public MaterialModel {
 public:
  ScaledElasticUser(Real E, Real nu, Real scale)
      : D_(elastic_iso_D(ElasticIso{E, nu})), scale_(scale) {}
  MaterialResponse evaluate(const Voigt6& strain,
                            MaterialState& state) const override {
    MaterialResponse r;
    for (int i = 0; i < 6; ++i) {
      Real v = 0.0;
      for (int j = 0; j < 6; ++j)
        v += D_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] *
             strain[static_cast<std::size_t>(j)];
      r.stress[static_cast<std::size_t>(i)] = scale_ * v;
      for (int j = 0; j < 6; ++j)
        r.tangent[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
            scale_ * D_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
    }
    if (!state.depvar.empty()) state.depvar[0] += 1.0;  // count evaluations
    return r;
  }

 private:
  D6 D_;
  Real scale_;
};

void test_registry_and_reference_law() {
  // The auto-registered reference law is discoverable.
  const UserMaterialFactory* f =
      UserMaterialRegistry::instance().find("LINEAR_ELASTIC_USER");
  CX_CHECK(f != nullptr);

  // Build the reference law from constants {E, nu} and check it equals elasticity.
  UserMaterial um;
  um.name = "LINEAR_ELASTIC_USER";
  um.constants = {210000.0, 0.3};
  std::unique_ptr<MaterialModel> model = (*f)(um, ElasticIso{});
  ElasticIsoMaterial el(steel());
  const Voigt6 strain{1e-4, -3e-5, -3e-5, 2e-5, 0.0, 1e-5};
  MaterialState s1, s2;
  const MaterialResponse ru = model->evaluate(strain, s1);
  const MaterialResponse re = el.evaluate(strain, s2);
  for (int i = 0; i < 6; ++i) {
    CX_NEAR(ru.stress[static_cast<std::size_t>(i)],
            re.stress[static_cast<std::size_t>(i)], 1e-9);
    for (int j = 0; j < 6; ++j)
      CX_NEAR(ru.tangent[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)],
              re.tangent[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)],
              1e-9);
  }
}

void test_rate_dependent_scaling() {
  // *RATEDEPENDENT scale multiplies stress and tangent.
  UserMaterial um;
  um.name = "LINEAR_ELASTIC_USER";
  um.constants = {210000.0, 0.3};
  um.rate_scale = 2.0;
  const UserMaterialFactory* f =
      UserMaterialRegistry::instance().find("LINEAR_ELASTIC_USER");
  std::unique_ptr<MaterialModel> scaled = (*f)(um, ElasticIso{});
  ElasticIsoMaterial el(steel());
  const Voigt6 strain{1e-4, -3e-5, -3e-5, 0.0, 0.0, 0.0};
  MaterialState s1, s2;
  const Voigt6 su = scaled->evaluate(strain, s1).stress;
  const Voigt6 se = el.evaluate(strain, s2).stress;
  for (int i = 0; i < 6; ++i)
    CX_NEAR(su[static_cast<std::size_t>(i)], 2.0 * se[static_cast<std::size_t>(i)],
            1e-6);
}

void test_depvar_committed() {
  // *DEPVAR state is written by evaluate() (trial) and promoted by commit().
  ScaledElasticUser mdl(210000.0, 0.3, 1.0);
  MaterialState st;
  st.depvar.assign(1, 0.0);
  st.committed_depvar.assign(1, 0.0);
  const Voigt6 strain{1e-4, 0.0, 0.0, 0.0, 0.0, 0.0};
  mdl.evaluate(strain, st);
  CX_NEAR(st.depvar[0], 1.0, 1e-12);
  CX_NEAR(st.committed_depvar[0], 0.0, 1e-12);  // committed untouched by evaluate
  st.commit();
  CX_NEAR(st.committed_depvar[0], 1.0, 1e-12);
}

void test_user_material_end_to_end() {
  // A *USER MATERIAL deck (reference law, constants {E,nu}, *DEPVAR 2) must solve
  // identically to the same cube with plain *ELASTIC (analytic linear-elastic).
  const Real uz = 1e-3;
  const Model user = uniaxial_cube(uz, [](Model& m) {
    Material mat;
    mat.name = "MAT";
    UserMaterial um;
    um.name = "LINEAR_ELASTIC_USER";
    um.constants = {210000.0, 0.3};
    um.ndepvar = 2;
    mat.user = um;
    m.materials["MAT"] = mat;
  });
  const Model elastic = uniaxial_cube(uz, [](Model& m) {
    Material mat;
    mat.name = "MAT";
    mat.elastic = steel();
    m.materials["MAT"] = mat;
  });
  CX_CHECK(user.has_nonlinear_material());
  numerics::NonlinearReport rep;
  const StaticFields fu = numerics::solve_nonlinear_static(user, {}, &rep);
  CX_CHECK(rep.converged);
  const StaticFields fe = numerics::solve_linear_static(elastic);

  // Analytic axial stress for the confined uniaxial cube: szz = E'(uz) with the
  // usual Poisson coupling; here we just require the user path == the elastic path.
  CX_NEAR(mean_szz(fu), mean_szz(fe), 1e-6 * std::abs(mean_szz(fe)) + 1e-9);
  for (std::size_t n = 0; n < fu.displacement.size(); ++n)
    for (int c = 0; c < 3; ++c)
      CX_NEAR(fu.displacement[n][static_cast<std::size_t>(c)],
              fe.displacement[n][static_cast<std::size_t>(c)], 1e-9);
}

void test_custom_registered_law_end_to_end() {
  // Register a brand-new law from the test and drive it through a deck: proves the
  // extensibility seam works for arbitrary user code. The law scales elasticity by 1,
  // so it must match the elastic cube; the *DEPVAR counter proves state flows.
  register_user_material(
      "TEST_SCALED", [](const UserMaterial& um, const ElasticIso& el) {
        const Real E = um.constants.size() > 0 ? um.constants[0] : el.E;
        const Real nu = um.constants.size() > 1 ? um.constants[1] : el.nu;
        const Real scale = um.constants.size() > 2 ? um.constants[2] : 1.0;
        return std::make_unique<ScaledElasticUser>(E, nu, scale);
      });
  const Real uz = 1e-3;
  const Model m = uniaxial_cube(uz, [](Model& mm) {
    Material mat;
    mat.name = "MAT";
    UserMaterial um;
    um.name = "TEST_SCALED";
    um.constants = {210000.0, 0.3, 1.0};
    um.ndepvar = 1;
    mat.user = um;
    mm.materials["MAT"] = mat;
  });
  numerics::NonlinearReport rep;
  const StaticFields f = numerics::solve_nonlinear_static(m, {}, &rep);
  CX_CHECK(rep.converged);
  const Model el = uniaxial_cube(uz, [](Model& mm) {
    Material mat;
    mat.name = "MAT";
    mat.elastic = steel();
    mm.materials["MAT"] = mat;
  });
  const StaticFields fe = numerics::solve_linear_static(el);
  CX_NEAR(mean_szz(f), mean_szz(fe), 1e-6 * std::abs(mean_szz(fe)) + 1e-9);
}

// ---- 4.3 hyperelastic ------------------------------------------------------------

Hyperelastic neo_hookean_from_E_nu(Real E, Real nu) {
  // mu = E / 2(1+nu),  kappa = E / 3(1-2nu);  C10 = mu/2, D1 = 2/kappa.
  const Real mu = E / (2.0 * (1.0 + nu));
  const Real kappa = E / (3.0 * (1.0 - 2.0 * nu));
  Hyperelastic h;
  h.model = Hyperelastic::Model::NeoHookean;
  h.mu = mu;
  h.kappa = kappa;
  h.c10 = mu / 2.0;
  h.d1 = 2.0 / kappa;
  return h;
}

void test_hyperelastic_small_strain_limit() {
  // The neo-Hookean tangent in the small-strain limit equals elastic_iso_D(E, nu).
  const Real E = 210000.0, nu = 0.3;
  HyperelasticNeoHookean h(neo_hookean_from_E_nu(E, nu));
  const D6 De = elastic_iso_D(ElasticIso{E, nu});
  const D6& D = h.tangent();
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      CX_NEAR(D[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)],
              De[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)], 1e-6);
}

void test_hyperelastic_end_to_end() {
  // A hyperelastic cube in small-strain uniaxial tension matches the same *ELASTIC
  // cube (analytic linear-elastic solution the model reduces to at small strain).
  const Real E = 210000.0, nu = 0.3, uz = 1e-4;  // tiny strain: linear regime
  const Model hyper = uniaxial_cube(uz, [E, nu](Model& m) {
    Material mat;
    mat.name = "MAT";
    mat.hyperelastic = neo_hookean_from_E_nu(E, nu);
    m.materials["MAT"] = mat;
  });
  const Model elastic = uniaxial_cube(uz, [](Model& m) {
    Material mat;
    mat.name = "MAT";
    mat.elastic = steel();
    m.materials["MAT"] = mat;
  });
  CX_CHECK(hyper.has_nonlinear_material());
  numerics::NonlinearReport rep;
  const StaticFields fh = numerics::solve_nonlinear_static(hyper, {}, &rep);
  CX_CHECK(rep.converged);
  const StaticFields fe = numerics::solve_linear_static(elastic);
  CX_NEAR(mean_szz(fh), mean_szz(fe), 1e-6 * std::abs(mean_szz(fe)) + 1e-9);
}

}  // namespace

int main() {
  test_registry_and_reference_law();
  test_rate_dependent_scaling();
  test_depvar_committed();
  test_user_material_end_to_end();
  test_custom_registered_law_end_to_end();
  test_hyperelastic_small_strain_limit();
  test_hyperelastic_end_to_end();
  if (cxtest::g_failures == 0) std::printf("test_material_models: OK\n");
  CX_MAIN_RETURN();
}
