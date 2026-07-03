#pragma once
#include <array>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "calculixpp/core/material.hpp"
#include "calculixpp/core/results.hpp"
#include "calculixpp/core/types.hpp"

// Material-point (constitutive) interface. A MaterialModel maps a strain increment
// (here: total strain, small-strain) plus a per-integration-point State to a Cauchy
// stress and a consistent tangent D_t, both in Voigt order xx,yy,zz,xy,xz,yz with
// engineering shear. Element kernels evaluate a MaterialModel at each Gauss point so
// nonlinear constitutive laws (plasticity, hyperelasticity — workstream 4) plug into
// the same B^T sigma / B^T D_t B assembly without touching the element code.
// (spec: material-models — material-point foundation.)
namespace cxpp::fem {

inline constexpr int kVoigt6 = 6;

// 6x6 tangent (constitutive) matrix, order xx,yy,zz,xy,xz,yz.
using D6 = std::array<std::array<Real, kVoigt6>, kVoigt6>;

// Per-integration-point history/state. Empty for linear elasticity; plasticity
// carries the plastic strain tensor, equivalent (accumulated) plastic strain, and
// the kinematic back stress. Passed by reference so evaluate() can advance the
// TRIAL history from the COMMITTED history at the current trial strain.
//
// Trial vs. committed: the Newton driver reassembles the tangent every iteration, so
// evaluate() must always return-map from the *converged* (committed) history, not
// from whatever the previous iteration wrote. It therefore reads `committed_*` and
// writes only the `*` (trial) fields; the driver calls commit() on the whole
// MaterialPoints store once an increment converges (spec: material-models — trial
// state per Newton iteration, committed on convergence).
struct MaterialState {
  // Committed (start-of-increment) history.
  Voigt6 committed_plastic_strain{};   // eps^p (engineering shear)
  Real committed_eqplastic{0.0};       // accumulated equivalent plastic strain
  Voigt6 committed_back_stress{};      // kinematic back stress alpha

  // Trial (current Newton iteration) history, written by evaluate().
  Voigt6 plastic_strain{};
  Real eqplastic{0.0};
  Voigt6 back_stress{};

  // Generic solution-dependent state variables for *USER MATERIAL / *DEPVAR
  // (spec: material-models 4.6). A user MaterialModel reads `committed_depvar`
  // (start-of-increment) and writes the advanced values into `depvar`; commit()
  // promotes them exactly like the plastic history. Sized by make_material_points
  // from the material's *DEPVAR count (empty for built-in models).
  std::vector<Real> committed_depvar;
  std::vector<Real> depvar;

  // Copy the trial history into the committed slot (called on increment convergence).
  void commit() {
    committed_plastic_strain = plastic_strain;
    committed_eqplastic = eqplastic;
    committed_back_stress = back_stress;
    committed_depvar = depvar;
  }
};

// Result of a constitutive evaluation at one integration point.
struct MaterialResponse {
  Voigt6 stress{};  // sigma
  D6 tangent{};     // D_t = d sigma / d strain
};

// Constitutive-law interface. `evaluate` is const (the model itself is stateless);
// all history lives in the mutable `state`.
class MaterialModel {
 public:
  virtual ~MaterialModel() = default;
  virtual MaterialResponse evaluate(const Voigt6& strain,
                                    MaterialState& state) const = 0;
};

// Linear isotropic elasticity: sigma = D * strain, tangent = D (constant), state
// unused. D is the same matrix elastic_iso_D produces, so a material-point assembly
// with this model reproduces the closed-form B^T D B element stiffness exactly.
class ElasticIsoMaterial final : public MaterialModel {
 public:
  explicit ElasticIsoMaterial(const ElasticIso& props);
  MaterialResponse evaluate(const Voigt6& strain,
                            MaterialState& state) const override;
  const D6& tangent() const { return D_; }

 private:
  D6 D_{};
};

// Rate-independent J2 (von Mises) plasticity with a radial-return mapping and an
// algorithmically CONSISTENT tangent (spec: material-models — return-mapping
// plasticity + consistent tangent). Small-strain, isotropic elasticity with
// piecewise-linear isotropic hardening (*PLASTIC), optional linear kinematic
// hardening (back stress), and *CYCLIC HARDENING. The elastic response and tangent
// reduce EXACTLY to ElasticIsoMaterial while the trial state stays inside the yield
// surface, so an unloaded/elastic model is unchanged.
//
// evaluate() reads the committed history from `state`, forms the elastic trial
// stress, and — if the trial exceeds the (isotropic + kinematic) yield surface —
// return-maps to the surface, writing the advanced plastic strain / eqplastic /
// back stress into the TRIAL slots of `state` (never the committed slots). It returns
// the updated Cauchy stress and the consistent tangent.
class J2PlasticMaterial final : public MaterialModel {
 public:
  J2PlasticMaterial(const ElasticIso& elastic, const Plastic& plastic);
  MaterialResponse evaluate(const Voigt6& strain,
                            MaterialState& state) const override;

 private:
  // Current yield stress and hardening slope at accumulated eqplastic `ep`, from the
  // piecewise-linear *PLASTIC table (flat extrapolation beyond the last point).
  void hardening(Real ep, Real& sigma_y, Real& h_iso) const;

  Real E_{};
  Real nu_{};
  Real mu_{};          // shear modulus
  Real h_kin_{};       // linear kinematic-hardening modulus
  std::vector<Real> ep_;      // plastic-strain abscissae (from *PLASTIC)
  std::vector<Real> sigma_y_; // yield ordinates
  D6 De_{};            // elastic tangent (order xx,yy,zz,xy,xz,yz)
};

// Compressible neo-Hookean hyperelasticity in the small-to-moderate-strain regime
// (spec: material-models 4.3 — *HYPERELASTIC). This is a *material-law* stand-in that
// plugs into the SAME small-strain B^T sigma / B^T D_t B assembly as the other models:
// it evaluates a decoupled deviatoric/volumetric response so a near-incompressible
// coefficient set (large kappa relative to mu, the u/p regime the spec flags) still
// gives a well-scaled tangent. It does NOT carry finite-strain kinematics (no F, no
// push-forward) — that is the deferred large-strain part of 4.3. For the small strains
// where it is valid the tangent reduces to the isotropic elastic D with Lame constants
// (lambda = kappa - 2/3 mu, mu), so a linear patch reproduces linear elasticity, and
// the model is validated against that analytic limit. State is unused (path-independent).
class HyperelasticNeoHookean final : public MaterialModel {
 public:
  explicit HyperelasticNeoHookean(const Hyperelastic& props);
  MaterialResponse evaluate(const Voigt6& strain,
                            MaterialState& state) const override;
  const D6& tangent() const { return D_; }

 private:
  Real mu_{};     // initial shear modulus
  Real kappa_{};  // initial bulk modulus
  D6 D_{};        // deviatoric+volumetric tangent (constant in the small-strain limit)
};

// ---------------------------------------------------------------------------------
// User-material plug-in registry (spec: material-models 4.6 — *USER MATERIAL C++
// interface). A *USER MATERIAL card in the deck names a factory registered here; the
// factory builds a MaterialModel from the parsed UserMaterial (constants, depvar
// count, rate scaling). This is the extensibility seam: a downstream user links a
// translation unit that calls register_user_material("MYLAW", ...) and their law is
// reachable from a deck without touching the solver core.
//
// The factory receives the fully-parsed UserMaterial plus the material's *ELASTIC
// block (may be empty) so simple laws can reuse the isotropic constants.
using UserMaterialFactory =
    std::function<std::unique_ptr<MaterialModel>(const UserMaterial&,
                                                 const ElasticIso&)>;

// Process-wide registry of user-material factories keyed by the *USER MATERIAL name
// (upper-cased). Thread-unsafe registration is expected at static-init / test setup.
class UserMaterialRegistry {
 public:
  static UserMaterialRegistry& instance();
  // Register (or replace) a factory under `name` (case-insensitive). Returns true.
  bool register_factory(const std::string& name, UserMaterialFactory factory);
  // Look up a factory; nullptr target if none registered under `name`.
  const UserMaterialFactory* find(const std::string& name) const;

 private:
  std::unordered_map<std::string, UserMaterialFactory> factories_;
};

// Convenience free function mirroring the class API (the documented registration hook).
bool register_user_material(const std::string& name, UserMaterialFactory factory);

// A built-in reference user material proving the plug-in path end to end
// (spec: material-models 4.6). It is a linear-elastic law driven ENTIRELY by the
// *USER MATERIAL constants: constants = {E, nu}. This demonstrates that a deck using
// *USER MATERIAL, CONSTANTS=2 / E,nu with *DEPVAR routes through the registry, builds
// this model, and reproduces the isotropic elastic result — a self-contained proof of
// the interface without shipping a third-party law. It also records a *DEPVAR slot
// (writes the equivalent strain into depvar[0] when present) and applies a
// *RATEDEPENDENT scale to the tangent+stress if one was given, exercising both fields.
class LinearUserMaterial final : public MaterialModel {
 public:
  LinearUserMaterial(Real E, Real nu, Real rate_scale);
  MaterialResponse evaluate(const Voigt6& strain,
                            MaterialState& state) const override;
  // Register this reference law under the given name (default "LINEAR_ELASTIC_USER").
  static bool register_default(
      const std::string& name = "LINEAR_ELASTIC_USER");

 private:
  D6 D_{};
  Real rate_scale_{1.0};
};

}  // namespace cxpp::fem
