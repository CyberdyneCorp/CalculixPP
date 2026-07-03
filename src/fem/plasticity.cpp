// Rate-independent J2 (von Mises) plasticity: radial-return mapping with a consistent
// (algorithmic) tangent. Math reference (reimplemented, not copied): Simo & Hughes,
// "Computational Inelasticity" (1998), Box 3.2 (return map) and Box 3.4 (consistent
// elastoplastic tangent); CalculiX umat_lin_iso_creep.f / plas_iso.f only as a
// behavioral oracle. Voigt order xx,yy,zz,xy,xz,yz with engineering shear.
//
// Internally the deviatoric algebra is done on the true (tensor) stress/strain so the
// von Mises norm and the return direction are correct despite the engineering-shear
// Voigt convention (the shear rows of the strain use gamma = 2*eps).
#include <array>
#include <cmath>
#include <stdexcept>

#include "calculixpp/fem/element.hpp"
#include "calculixpp/fem/material_model.hpp"

namespace cxpp::fem {
namespace {

// sqrt(2/3) and sqrt(3/2): the J2 factors linking the deviatoric tensor norm ||s||
// to the equivalent (von Mises) stress q = sqrt(3/2) ||s|| and the equivalent plastic
// strain increment dep = sqrt(2/3) ||de^p||.
constexpr Real kSqrt2_3 = 0.816496580927726;  // sqrt(2/3)
constexpr Real kSqrt3_2 = 1.224744871391589;  // sqrt(3/2)

// Double contraction s:s for a symmetric tensor stored in Voigt (xx,yy,zz,xy,xz,yz).
// The off-diagonal (shear) components appear twice in the full 3x3 tensor, so they
// are weighted by 2. Valid for a STRESS-like Voigt vector (true shear components).
Real contract(const Voigt6& s) {
  return s[0] * s[0] + s[1] * s[1] + s[2] * s[2] +
         2.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]);
}

// Deviatoric part of a stress-like Voigt vector (removes the mean of the 3 normals).
Voigt6 deviator(const Voigt6& t) {
  const Real p = (t[0] + t[1] + t[2]) / 3.0;
  return {t[0] - p, t[1] - p, t[2] - p, t[3], t[4], t[5]};
}

}  // namespace

J2PlasticMaterial::J2PlasticMaterial(const ElasticIso& elastic,
                                     const Plastic& plastic)
    : E_(elastic.E), nu_(elastic.nu) {
  if (plastic.empty())
    throw std::runtime_error("J2PlasticMaterial needs a non-empty *PLASTIC table");
  mu_ = E_ / (2.0 * (1.0 + nu_));
  h_kin_ = plastic.kinematic_modulus;
  ep_ = plastic.eqplastic;
  sigma_y_ = plastic.yield;
  // Kinematic / combined hardening: derive the (linear) back-stress modulus H_kin
  // from the table's first slope when the parser did not supply it directly. For
  // pure KINEMATIC hardening the yield surface does not grow (isotropic yield stays
  // at the initial value); the table slope drives the back stress instead, so we
  // flatten the isotropic curve to the initial yield (hardening() then returns a
  // constant sigma_y with zero isotropic slope).
  if (plastic.hardening == Plastic::Hardening::Kinematic ||
      plastic.hardening == Plastic::Hardening::Combined) {
    if (h_kin_ == 0.0 && ep_.size() >= 2 && ep_[1] > ep_[0])
      h_kin_ = (sigma_y_[1] - sigma_y_[0]) / (ep_[1] - ep_[0]);
    if (plastic.hardening == Plastic::Hardening::Kinematic)
      for (Real& y : sigma_y_) y = sigma_y_.front();
  }
  De_ = elastic_iso_D(elastic);
}

void J2PlasticMaterial::hardening(Real ep, Real& sigma_y, Real& h_iso) const {
  const std::size_t n = ep_.size();
  if (n == 1 || ep <= ep_.front()) {  // perfect plasticity or at/below the first pt
    sigma_y = sigma_y_.front();
    h_iso = 0.0;
    return;
  }
  if (ep >= ep_.back()) {  // flat extrapolation beyond the last table point
    sigma_y = sigma_y_.back();
    h_iso = 0.0;
    return;
  }
  for (std::size_t i = 0; i + 1 < n; ++i) {
    if (ep >= ep_[i] && ep <= ep_[i + 1]) {
      h_iso = (sigma_y_[i + 1] - sigma_y_[i]) / (ep_[i + 1] - ep_[i]);
      sigma_y = sigma_y_[i] + h_iso * (ep - ep_[i]);
      return;
    }
  }
  sigma_y = sigma_y_.back();
  h_iso = 0.0;
}

MaterialResponse J2PlasticMaterial::evaluate(const Voigt6& strain,
                                             MaterialState& state) const {
  // Elastic trial: stress from the elastic strain (total minus committed plastic).
  // Note strain uses engineering shear; converting the plastic strain back to the
  // same convention keeps the elastic law sigma = De*(eps - eps^p) consistent.
  Voigt6 el_strain{};
  for (int i = 0; i < kVoigt6; ++i)
    el_strain[static_cast<std::size_t>(i)] =
        strain[static_cast<std::size_t>(i)] -
        state.committed_plastic_strain[static_cast<std::size_t>(i)];
  Voigt6 sig_tr{};
  for (int i = 0; i < kVoigt6; ++i) {
    Real v = 0.0;
    for (int j = 0; j < kVoigt6; ++j)
      v += De_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] *
           el_strain[static_cast<std::size_t>(j)];
    sig_tr[static_cast<std::size_t>(i)] = v;
  }

  // Relative deviatoric trial stress xi = dev(sigma_tr) - back_stress.
  const Voigt6 s_tr = deviator(sig_tr);
  Voigt6 xi{};
  for (int i = 0; i < kVoigt6; ++i)
    xi[static_cast<std::size_t>(i)] =
        s_tr[static_cast<std::size_t>(i)] -
        state.committed_back_stress[static_cast<std::size_t>(i)];
  const Real xi_norm = std::sqrt(contract(xi));   // ||xi||
  const Real q_tr = kSqrt3_2 * xi_norm;           // von Mises trial stress

  Real sy0 = 0.0, h0 = 0.0;
  hardening(state.committed_eqplastic, sy0, h0);
  const Real f_tr = q_tr - sy0;  // yield function at the trial state

  MaterialResponse r;
  if (f_tr <= 0.0 || xi_norm == 0.0) {  // ELASTIC step: nothing plastic advances
    r.stress = sig_tr;
    r.tangent = De_;
    state.plastic_strain = state.committed_plastic_strain;
    state.eqplastic = state.committed_eqplastic;
    state.back_stress = state.committed_back_stress;
    return r;
  }

  // PLASTIC step: solve q_tr - 3 mu dgamma - sigma_y(ep0+dgamma) - H_kin dgamma = 0
  // for the plastic multiplier dgamma (>= 0) by a local Newton iteration. The
  // piecewise-linear isotropic hardening makes h_iso depend on ep = ep0 + dgamma.
  const Real ep0 = state.committed_eqplastic;
  Real dgamma = 0.0;
  for (int it = 0; it < 50; ++it) {
    Real sy = 0.0, h_iso = 0.0;
    hardening(ep0 + dgamma, sy, h_iso);
    const Real g = q_tr - 3.0 * mu_ * dgamma - h_kin_ * dgamma - sy;
    const Real dg = -(3.0 * mu_ + h_kin_ + h_iso);  // dg/d(dgamma)
    const Real step = -g / dg;
    dgamma += step;
    if (dgamma < 0.0) dgamma = 0.0;
    if (std::abs(step) <= 1e-14 * (1.0 + std::abs(dgamma))) break;
  }

  // Return direction n = xi / ||xi|| (unit deviatoric tensor). The plastic strain
  // increment (tensor) is dgamma * sqrt(3/2) * n; in engineering-shear Voigt the
  // shear rows double.
  Voigt6 n{};
  for (int i = 0; i < kVoigt6; ++i)
    n[static_cast<std::size_t>(i)] = xi[static_cast<std::size_t>(i)] / xi_norm;

  // Updated deviatoric stress s = s_tr - (2 mu dgamma sqrt(3/2)) n. The corrector
  // magnitude 2 mu dgamma sqrt(3/2) along n makes q = sqrt(3/2)||dev|| drop by exactly
  // 3 mu dgamma, matching the return-map residual solved above.
  const Real dev_scale = 2.0 * mu_ * dgamma * kSqrt3_2;
  // Back-stress increment magnitude along n: dalpha = (2/3) H_kin de^p_tensor, with
  // de^p_tensor = dgamma sqrt(3/2) n  ->  alpha_scale = (2/3) H_kin dgamma sqrt(3/2).
  const Real alpha_scale = (2.0 / 3.0) * h_kin_ * dgamma * kSqrt3_2;
  const Real mean = (sig_tr[0] + sig_tr[1] + sig_tr[2]) / 3.0;
  Voigt6 sigma{};
  for (int i = 0; i < kVoigt6; ++i) {
    const Real dev = s_tr[static_cast<std::size_t>(i)] -
                     dev_scale * n[static_cast<std::size_t>(i)];
    sigma[static_cast<std::size_t>(i)] =
        dev + (i < 3 ? mean : 0.0);  // add back the (unchanged) hydrostatic part
  }
  r.stress = sigma;

  // Advance the TRIAL history (committed only by the driver on convergence).
  // Plastic strain increment tensor de^p = dgamma*sqrt(3/2)*n; engineering shear
  // doubles the shear rows.
  const Real dep_scale = dgamma * kSqrt3_2;
  for (int i = 0; i < kVoigt6; ++i) {
    const Real dep_i =
        dep_scale * n[static_cast<std::size_t>(i)] * (i < 3 ? 1.0 : 2.0);
    state.plastic_strain[static_cast<std::size_t>(i)] =
        state.committed_plastic_strain[static_cast<std::size_t>(i)] + dep_i;
    state.back_stress[static_cast<std::size_t>(i)] =
        state.committed_back_stress[static_cast<std::size_t>(i)] +
        alpha_scale * n[static_cast<std::size_t>(i)];
  }
  state.eqplastic = ep0 + dgamma;

  // Consistent (algorithmic) tangent (Simo & Hughes, Box 3.4). Written directly in
  // the engineering-shear stress/strain Voigt convention used everywhere here.
  //
  // The elastic tangent De maps engineering-shear strain to stress. The plastic
  // corrector only touches the deviatoric response, so we START from De and subtract
  // two rank-structured terms:
  //
  //   D_alg = De - beta_dev * P - beta_n * (n n^T)
  //
  // where n is the (stress-like) unit deviatoric direction and P is the engineering
  // deviatoric projector applied on the shear-halving of the strain columns. The key
  // convention fact: for a stress-like deviator n and an engineering strain vector
  // de, the plain Voigt dot product n^T de EQUALS the tensor contraction n:de_tensor
  // (the shear doubling of engineering strain cancels the tensor factor 2). So the
  // outer product n n^T needs no extra shear factor.
  //
  // Scalars (Simo-Hughes):
  //   ||xi||        = deviatoric relative-stress norm (tensor)
  //   theta    = 1 - 2 mu dgamma sqrt(3/2) / ||xi||        (deviatoric softening)
  //   theta_bar= 1/(1 + (H_iso + H_kin)/(3 mu)) - (1 - theta)
  //   D_alg = De - 2 mu (1 - theta) * P_dev_eng - 2 mu theta_bar * (n n^T)
  Real sy = 0.0, h_iso = 0.0;
  hardening(ep0 + dgamma, sy, h_iso);
  const Real theta = 1.0 - (2.0 * mu_ * dgamma * kSqrt3_2) / xi_norm;
  const Real denom = 1.0 + (h_iso + h_kin_) / (3.0 * mu_);
  const Real theta_bar = 1.0 / denom - (1.0 - theta);

  // Idev = 2 mu * I_dev in engineering-shear Voigt: normals get 2mu*(delta_ij - 1/3),
  // shears get mu (matching De's shear diagonal). Subtracting 2 mu (1-theta) I_dev and
  // 2 mu theta_bar (n n^T) from De yields the algorithmic tangent.
  D6 Idev{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      Idev[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
          2.0 * mu_ * ((i == j ? 1.0 : 0.0) - 1.0 / 3.0);
  for (int i = 3; i < 6; ++i)
    Idev[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] = mu_;

  D6 D{};
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 6; ++j)
      D[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
          De_[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] -
          (1.0 - theta) * Idev[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] -
          2.0 * mu_ * theta_bar * n[static_cast<std::size_t>(i)] *
              n[static_cast<std::size_t>(j)];
  r.tangent = D;
  return r;
}

}  // namespace cxpp::fem
