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

// A temperature-dependent scalar material property, stored as a piecewise-linear
// table of (value, temperature) rows in strictly increasing temperature (spec:
// material-models — temperature-dependent *CONDUCTIVITY / *SPECIFIC HEAT / *EXPANSION).
// A CalculiX property card lists one row per data line, each ending in a temperature;
// a single (constant) row is the common case. `at(T)` interpolates piecewise-linearly
// and CLAMPS below the first / above the last tabulated temperature (CalculiX's
// convention). A single-row table returns that constant for every T, so a constant
// property is byte-for-byte the pre-table behavior.
struct PropertyTable {
  std::vector<Real> value;  // property ordinates (k, c, or alpha), one per row
  std::vector<Real> temp;   // temperature abscissae, strictly increasing

  PropertyTable() = default;
  // Constant (single-row) table: the common case, and a convenient implicit
  // conversion so `Thermal{k, c}` / `Expansion{alpha, tref}` keep working.
  PropertyTable(Real constant) : value{constant}, temp{0.0} {}

  bool empty() const { return value.empty(); }
  Real first() const { return value.empty() ? 0.0 : value.front(); }

  // Piecewise-linear interpolation at temperature T, clamped outside the table.
  Real at(Real T) const {
    if (value.empty()) return 0.0;
    if (value.size() == 1 || T <= temp.front()) return value.front();
    if (T >= temp.back()) return value.back();
    std::size_t i = 1;
    while (i < temp.size() && T > temp[i]) ++i;
    const Real t0 = temp[i - 1], t1 = temp[i];
    const Real w = (t1 == t0) ? 0.0 : (T - t0) / (t1 - t0);
    return value[i - 1] + w * (value[i] - value[i - 1]);
  }
};

// Thermal material properties for heat-transfer analysis (spec:
// heat-transfer-analysis / material-models — *CONDUCTIVITY, *SPECIFIC HEAT).
// Conductivity is stored isotropically as a temperature-dependent table used in the
// conduction element matrix Kt_e = ∫ ∇N_a·(k(T) ∇N_b) dV. Specific heat and density
// (via Material::density) feed the transient capacitance C_e = ∫ rho*c(T) N_a N_b dV.
// A single-row table is a constant property (byte-for-byte the pre-table path).
struct Thermal {
  PropertyTable conductivity;   // *CONDUCTIVITY (isotropic scalar k(T))
  PropertyTable specific_heat;  // *SPECIFIC HEAT (c(T)), used by transient capacitance
  bool empty() const { return conductivity.empty() && specific_heat.empty(); }
};

// Thermal-expansion data (spec: heat-transfer-analysis — coupled; material-models —
// *EXPANSION). Isotropic linear expansion: the thermal strain is
// eps_th = alpha (T - Tref) applied to the three normal components (Voigt xx,yy,zz),
// which the mechanical path subtracts from the total strain so
// sigma = D (eps_mech - eps_th). `alpha` is the *EXPANSION coefficient (isotropic,
// temperature-dependent table evaluated at the point temperature); `t_ref` is the
// ZERO parameter on *EXPANSION (the reference / stress-free temperature, default 0).
// A single-row table is a constant coefficient (byte-for-byte the pre-table path).
struct Expansion {
  PropertyTable alpha;  // *EXPANSION coefficient (isotropic alpha(T))
  Real t_ref{0.0};      // reference (stress-free) temperature — *EXPANSION, ZERO=
  bool empty() const { return alpha.empty(); }
};

struct Material {
  std::string name;
  std::optional<ElasticIso> elastic;
  std::optional<Real> density;
  std::optional<Plastic> plastic;  // *PLASTIC / *CYCLIC HARDENING (J2 plasticity)
  std::optional<Hyperelastic> hyperelastic;  // *HYPERELASTIC / *HYPERFOAM (4.3)
  std::optional<UserMaterial> user;          // *USER MATERIAL / *DEPVAR (4.6)
  std::optional<Thermal> thermal;            // *CONDUCTIVITY / *SPECIFIC HEAT (Phase 3)
  std::optional<Expansion> expansion;        // *EXPANSION (thermal-strain coupling)

  // Effective isotropic elastic properties for the linear DOF-map / initial-tangent
  // assembly. Returns the *ELASTIC block if present; otherwise derives an equivalent
  // (E, nu) from a *HYPERELASTIC neo-Hookean (mu, kappa) or *USER MATERIAL constants
  // {E, nu}, so a purely nonlinear material still yields a valid linear seed. Returns
  // nullopt only when the material carries none of these. Defined in model.cpp.
  std::optional<ElasticIso> effective_elastic() const;
};

}  // namespace cxpp
