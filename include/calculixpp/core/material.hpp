#pragma once
#include <optional>
#include <string>
#include <vector>

#include "calculixpp/core/types.hpp"

namespace cxpp {

// Isotropic linear elasticity (spec: material-models — Phase 1).
struct ElasticIso {
  Real E{};   // Young's modulus
  Real nu{};  // Poisson's ratio
};

// Rate-independent J2 (von Mises) plasticity data (spec: material-models —
// "Return-mapping rate-independent plasticity: *PLASTIC isotropic + kinematic,
// *CYCLIC HARDENING"). The hardening curve is a piecewise-linear table of yield
// stress vs. equivalent plastic strain (*PLASTIC), optionally combined with linear
// kinematic hardening (back stress) and *CYCLIC HARDENING data.
//
// Voigt order elsewhere is xx,yy,zz,xy,xz,yz (engineering shear). Plasticity keeps
// its own accumulated plastic strain per integration point (see MaterialState).
struct Plastic {
  // Selects how the total hardening is split between isotropic (grows the yield
  // surface) and kinematic (translates it via a back stress). CalculiX *PLASTIC
  // HARDENING=ISOTROPIC (default) / KINEMATIC / COMBINED.
  enum class Hardening { Isotropic, Kinematic, Combined };
  Hardening hardening{Hardening::Isotropic};

  // Piecewise-linear isotropic hardening curve: yield[i] at eqplastic[i], strictly
  // increasing plastic strain, yield[0] is the initial yield stress at zero plastic
  // strain. A single point means perfect plasticity (constant yield). For KINEMATIC
  // hardening this table (via *CYCLIC HARDENING or *PLASTIC) supplies the back-stress
  // modulus through its slope.
  std::vector<Real> eqplastic;  // equivalent plastic strain abscissae
  std::vector<Real> yield;      // yield stress ordinates

  // Linear kinematic-hardening modulus H_kin (back stress alpha grows as
  // (2/3) H_kin * plastic strain). Nonzero only for KINEMATIC / COMBINED. When zero,
  // the model is purely isotropic. Derived from *CYCLIC HARDENING / *PLASTIC slope.
  Real kinematic_modulus{0.0};

  bool empty() const { return yield.empty(); }
};

// Hyperelastic / foam material data (spec: material-models 4.3 — *HYPERELASTIC,
// *HYPERFOAM). Small-to-moderate-strain compressible neo-Hookean is the model
// implemented here; the coefficient convention follows CalculiX's *HYPERELASTIC,
// N=1 (C10, D1): strain-energy W = C10 (I1bar - 3) + (1/D1)(J - 1)^2, with the
// initial shear modulus mu0 = 2*C10 and initial bulk modulus K0 = 2/D1. The
// `mu`/`kappa` fields hold those derived initial moduli (large-strain finite
// kinematics deferred — see tasks.md 4.3 note).
struct Hyperelastic {
  enum class Model { NeoHookean };
  Model model{Model::NeoHookean};
  Real c10{0.0};  // *HYPERELASTIC neo-Hookean C10
  Real d1{0.0};   // *HYPERELASTIC neo-Hookean D1 (compressibility)
  Real mu{0.0};   // derived initial shear modulus  = 2*C10
  Real kappa{0.0};  // derived initial bulk modulus = 2/D1 (0 -> incompressible)
  bool empty() const { return mu == 0.0 && c10 == 0.0; }
};

// User-defined material (spec: material-models 4.6 — *USER MATERIAL C++ interface,
// *DEPVAR state, *RATEDEPENDENT scaling). The deck names a `name` that must match a
// MaterialModel factory registered in the user-material registry (see
// fem/material_model.hpp UserMaterialRegistry). `constants` are the *USER MATERIAL
// data-line values (CONSTANTS=n on the card), passed verbatim to the factory.
// `ndepvar` is the *DEPVAR state-variable count (per integration point). A nonzero
// `rate_scale` records a *RATEDEPENDENT scaling factor the model may apply.
struct UserMaterial {
  std::string name;                 // registry key (*USER MATERIAL, [NAME=]... )
  std::vector<Real> constants;      // *USER MATERIAL data values
  int ndepvar{0};                   // *DEPVAR state-variable count
  std::optional<Real> rate_scale;   // *RATEDEPENDENT scaling factor (if given)
  bool empty() const { return name.empty(); }
};

struct Material {
  std::string name;
  std::optional<ElasticIso> elastic;
  std::optional<Real> density;
  std::optional<Plastic> plastic;  // *PLASTIC / *CYCLIC HARDENING (J2 plasticity)
  std::optional<Hyperelastic> hyperelastic;  // *HYPERELASTIC / *HYPERFOAM (4.3)
  std::optional<UserMaterial> user;          // *USER MATERIAL / *DEPVAR (4.6)

  // Effective isotropic elastic properties for the linear DOF-map / initial-tangent
  // assembly. Returns the *ELASTIC block if present; otherwise derives an equivalent
  // (E, nu) from a *HYPERELASTIC neo-Hookean (mu, kappa) or *USER MATERIAL constants
  // {E, nu}, so a purely nonlinear material still yields a valid linear seed. Returns
  // nullopt only when the material carries none of these. Defined in model.cpp.
  std::optional<ElasticIso> effective_elastic() const;
};

}  // namespace cxpp
