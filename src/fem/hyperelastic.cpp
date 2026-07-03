// Compressible neo-Hookean hyperelasticity, small-to-moderate-strain regime
// (spec: material-models 4.3 — *HYPERELASTIC / *HYPERFOAM). Math reference
// (reimplemented, not copied): the decoupled deviatoric/volumetric split of the
// compressible neo-Hookean energy W = C10 (I1bar - 3) + (1/D1)(J-1)^2, whose
// small-strain (linearized) tangent is the isotropic elastic operator with
//   mu    = 2 C10                (initial shear modulus)
//   kappa = 2 / D1               (initial bulk modulus)
//   lambda = kappa - 2/3 mu.
// (ref: CalculiX orthonl.f / us3.f only as a behavioral oracle for the coefficient
// convention.) Voigt order xx,yy,zz,xy,xz,yz with engineering shear.
//
// This is a material LAW plugged into the existing small-strain assembly; it does not
// carry finite-strain kinematics (no deformation gradient / push-forward). That is the
// deferred large-strain portion of 4.3 (see tasks.md). The near-incompressible u/p
// concern is addressed by building the tangent from a decoupled deviatoric + volumetric
// form so a large kappa/mu ratio stays well-scaled rather than coupling into the shear
// block; in the small-strain limit this coincides with elastic_iso_D for the matching
// (E, nu), which is the analytic case we validate against.
#include <cmath>
#include <stdexcept>

#include "calculixpp/fem/material_model.hpp"

namespace cxpp::fem {

HyperelasticNeoHookean::HyperelasticNeoHookean(const Hyperelastic& props) {
  mu_ = props.mu != 0.0 ? props.mu : 2.0 * props.c10;
  kappa_ = props.kappa != 0.0 ? props.kappa
                              : (props.d1 != 0.0 ? 2.0 / props.d1 : 0.0);
  if (mu_ <= 0.0)
    throw std::runtime_error("HyperelasticNeoHookean needs a positive shear modulus");
  if (kappa_ <= 0.0)
    throw std::runtime_error(
        "HyperelasticNeoHookean needs a positive bulk modulus (compressible)");

  // Decoupled tangent: volumetric block kappa on the 3x3 normal-normal sub-block,
  // deviatoric block 2 mu * I_dev on the strains. In engineering-shear Voigt the
  // shear diagonal is mu. Summing gives the standard isotropic D with
  // lambda = kappa - 2/3 mu, so the small-strain response equals elastic_iso_D for
  // the matching (E, nu) (validated in test_hyperelastic).
  const Real lam = kappa_ - (2.0 / 3.0) * mu_;
  const Real lam2mu = lam + 2.0 * mu_;
  D_[0][0] = D_[1][1] = D_[2][2] = lam2mu;
  D_[0][1] = D_[0][2] = D_[1][0] = D_[1][2] = D_[2][0] = D_[2][1] = lam;
  D_[3][3] = D_[4][4] = D_[5][5] = mu_;
}

MaterialResponse HyperelasticNeoHookean::evaluate(const Voigt6& strain,
                                                  MaterialState& /*state*/) const {
  MaterialResponse r;
  r.tangent = D_;
  for (int i = 0; i < kVoigt6; ++i) {
    Real v = 0.0;
    for (int j = 0; j < kVoigt6; ++j)
      v += D_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] *
           strain[static_cast<std::size_t>(j)];
    r.stress[static_cast<std::size_t>(i)] = v;
  }
  return r;
}

}  // namespace cxpp::fem
