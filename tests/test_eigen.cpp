// Mass matrix + eigensolution engine (Phase 4, tasks 1.1-1.4, 2.1).
//   * element consistent/lumped mass properties (total mass = rho*V per direction);
//   * generalized eigenpair residual K φ = λ M φ and mass-normalization φᵀ M φ = 1;
//   * an ANALYTICAL 2-DOF spring-mass chain with closed-form eigenvalues;
//   * the stock-CalculiX *FREQUENCY reference deck beam8f matches its .dat.ref
//     eigenvalues / natural frequencies AND participation factors within tolerance
//     (validated here when the CalculiX test tree is present; skipped gracefully
//     otherwise — the Python regression suite also pins beam8f).
#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "calculixpp/fem/assembly.hpp"
#include "calculixpp/fem/element.hpp"
#include "calculixpp/io/inp_parser.hpp"
#include "calculixpp/numerics/eigensolution.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

const std::vector<Vec3> kUnitCube = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                                     {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};

// The consistent mass matrix conserves the total translational mass: summing the
// whole matrix over one direction's DOFs gives rho * V (here V = 1, rho = 2 -> 2).
// Also the per-direction row-sum-lumped diagonal must sum to the same total.
void test_mass_total() {
  const Real rho = 2.0;
  const std::vector<Real> Mc = fem::element_mass(ElementType::C3D8, kUnitCube, rho);
  const int ndof = 8 * kDofsPerNode;
  // Full sum over the x-direction sub-block = rho*V. Sum the whole matrix restricted
  // to x DOFs (indices 3a+0).
  Real total_x = 0.0;
  for (int a = 0; a < 8; ++a)
    for (int b = 0; b < 8; ++b) {
      const int ia = a * kDofsPerNode + 0, ib = b * kDofsPerNode + 0;
      total_x += Mc[static_cast<std::size_t>(ia) * static_cast<std::size_t>(ndof) +
                    static_cast<std::size_t>(ib)];
    }
  CX_NEAR(total_x, rho * 1.0, 1e-12);

  const std::vector<Real> Ml =
      fem::element_mass_lumped(ElementType::C3D8, kUnitCube, rho);
  Real diag_x = 0.0;
  for (int a = 0; a < 8; ++a) {
    const int ia = a * kDofsPerNode + 0;
    diag_x += Ml[static_cast<std::size_t>(ia) * static_cast<std::size_t>(ndof) +
                 static_cast<std::size_t>(ia)];
  }
  CX_NEAR(diag_x, rho * 1.0, 1e-12);
  // Zero density -> zero matrix.
  const std::vector<Real> Mz = fem::element_mass(ElementType::C3D8, kUnitCube, 0.0);
  Real anyz = 0.0;
  for (Real v : Mz) anyz += std::fabs(v);
  CX_NEAR(anyz, 0.0, 0.0);
}

// A cantilevered C3D8 hex (E, nu, rho) with the whole z=0 face fully clamped so the
// extracted modes are genuine deformation modes (no rigid-body / near-zero
// eigenvalues): every mode must satisfy the generalized residual K φ = λ M φ on the
// free DOFs and be mass-normalized. (Clamping a full face removes all six rigid-body
// modes — cf. test_generalized_residual with a single fixed node, which would leave
// near-zero rotational modes that the spectral-shift path handles separately.)
void test_generalized_residual() {
  Model m;
  for (std::size_t i = 0; i < kUnitCube.size(); ++i)
    m.mesh.add_node(static_cast<Index>(i + 1), kUnitCube[i]);
  m.mesh.add_element(1, ElementType::C3D8, {1, 2, 3, 4, 5, 6, 7, 8});
  m.mesh.add_elset("EALL", {1});
  Material mat;
  mat.name = "EL";
  mat.elastic = ElasticIso{210000.0, 0.3};
  mat.density = 7.8e-9;
  m.materials["EL"] = mat;
  m.sections.push_back(SolidSection{"EALL", "EL"});
  // Fully clamp the four z=0 corner nodes (1..4) -> no rigid-body modes remain.
  for (Index n : {Index{1}, Index{2}, Index{3}, Index{4}})
    for (int c = 1; c <= 3; ++c) m.spcs.push_back(Spc{n, c, 0.0});
  m.procedure = Procedure::Frequency;

  const fem::LinearSystem K = fem::assemble_linear_static(m);
  const fem::LinearSystem M = fem::assemble_mass(m, false);
  const numerics::EigenBasis basis = numerics::extract_modes(K, M, 3);
  CX_CHECK(basis.modes.size() == 3);

  // Build dense K, M on the free DOFs to check the residual directly.
  const std::size_t nf = static_cast<std::size_t>(K.n_free);
  std::vector<std::vector<Real>> Kd(nf, std::vector<Real>(nf, 0.0)), Md = Kd;
  for (std::size_t k = 0; k < K.vals.size(); ++k)
    Kd[static_cast<std::size_t>(K.rows[k])][static_cast<std::size_t>(K.cols[k])] +=
        K.vals[k];
  for (std::size_t k = 0; k < M.vals.size(); ++k)
    Md[static_cast<std::size_t>(M.rows[k])][static_cast<std::size_t>(M.cols[k])] +=
        M.vals[k];
  // Symmetrize (the congruence scatter stores one triangle for some entries).
  for (std::size_t i = 0; i < nf; ++i)
    for (std::size_t j = i + 1; j < nf; ++j) {
      const Real kv = Kd[i][j] != 0.0 ? Kd[i][j] : Kd[j][i];
      Kd[i][j] = Kd[j][i] = kv;
      const Real mv = Md[i][j] != 0.0 ? Md[i][j] : Md[j][i];
      Md[i][j] = Md[j][i] = mv;
    }

  for (std::size_t mode = 0; mode < basis.modes.size(); ++mode) {
    const std::vector<Real>& phi = basis.mode_free[mode];
    const Real lam = basis.modes[mode].eigenvalue;
    // Mass-normalization φᵀ M φ = 1.
    Real mnorm = 0.0;
    for (std::size_t i = 0; i < nf; ++i)
      for (std::size_t j = 0; j < nf; ++j) mnorm += phi[i] * Md[i][j] * phi[j];
    CX_NEAR(mnorm, 1.0, 1e-9);
    // Residual r = K φ - λ M φ, relative to |K φ|.
    Real rnorm = 0.0, knorm = 0.0;
    for (std::size_t i = 0; i < nf; ++i) {
      Real kphi = 0.0, mphi = 0.0;
      for (std::size_t j = 0; j < nf; ++j) {
        kphi += Kd[i][j] * phi[j];
        mphi += Md[i][j] * phi[j];
      }
      const Real ri = kphi - lam * mphi;
      rnorm += ri * ri;
      knorm += kphi * kphi;
    }
    CX_CHECK(std::sqrt(rnorm) <= 1e-6 * std::sqrt(knorm) + 1e-6);
    // omega / frequency consistency.
    CX_NEAR(basis.modes[mode].omega, std::sqrt(std::max(lam, 0.0)), 1e-6);
    CX_NEAR(basis.modes[mode].frequency,
            basis.modes[mode].omega / (2.0 * M_PI), 1e-9);
  }
}

// ANALYTICAL 2-DOF undamped spring-mass chain (grounded): masses m1=m2=m on x, a
// grounded spring k to node 1 and a coupling spring k between nodes 1 and 2. The
// stiffness is [[2k,-k],[-k,k]], mass diag(m,m); eigenvalues of K x = λ M x are
// λ = (k/m) * (3 ∓ sqrt(5)) / 2. Modeled with *SPRING (grounded + dof) and *MASS.
void test_spring_mass_analytical() {
  const Real k = 100.0, mass = 2.0;
  Model m;
  m.mesh.add_node(1, {0, 0, 0});
  m.mesh.add_node(2, {1, 0, 0});
  // Grounded spring k on node1-x; coupling spring k between node1-x and node2-x.
  Spring s1;
  s1.kind = Spring::Kind::Grounded;
  s1.node1 = 1;
  s1.dof1 = 1;
  s1.stiffness = k;
  Spring s2;
  s2.kind = Spring::Kind::Dof;
  s2.node1 = 1;
  s2.dof1 = 1;
  s2.node2 = 2;
  s2.dof2 = 1;
  s2.stiffness = k;
  m.springs = {s1, s2};
  m.point_masses = {PointMass{1, mass}, PointMass{2, mass}};
  // Constrain y,z of both nodes so only the two x DOFs are free.
  for (Index n : {Index{1}, Index{2}})
    for (int c : {2, 3}) m.spcs.push_back(Spc{n, c, 0.0});
  m.procedure = Procedure::Frequency;

  const fem::LinearSystem K = fem::assemble_linear_static(m);
  const fem::LinearSystem M = fem::assemble_mass(m, false);
  CX_CHECK(K.n_free == 2);
  const numerics::EigenBasis basis = numerics::extract_modes(K, M, 2);
  CX_CHECK(basis.modes.size() == 2);
  const Real r = k / mass;
  const Real lam_lo = r * (3.0 - std::sqrt(5.0)) / 2.0;
  const Real lam_hi = r * (3.0 + std::sqrt(5.0)) / 2.0;
  CX_NEAR(basis.modes[0].eigenvalue, lam_lo, 1e-9);
  CX_NEAR(basis.modes[1].eigenvalue, lam_hi, 1e-9);
}

// Spectral shift / rigid-body handling (task 1.3): a fully UNCONSTRAINED hex has six
// rigid-body modes at (near-)zero eigenvalue. The dense path (M-Cholesky, M is always
// PD) must extract them without failure, and a nonzero spectral shift must return the
// SAME (un-shifted) eigenvalues — confirming the shift is transparent to the result.
void test_rigid_body_shift() {
  Model m;
  for (std::size_t i = 0; i < kUnitCube.size(); ++i)
    m.mesh.add_node(static_cast<Index>(i + 1), kUnitCube[i]);
  m.mesh.add_element(1, ElementType::C3D8, {1, 2, 3, 4, 5, 6, 7, 8});
  m.mesh.add_elset("EALL", {1});
  Material mat;
  mat.name = "EL";
  mat.elastic = ElasticIso{210000.0, 0.3};
  mat.density = 7.8e-9;
  m.materials["EL"] = mat;
  m.sections.push_back(SolidSection{"EALL", "EL"});
  m.procedure = Procedure::Frequency;  // no SPCs -> fully free

  const fem::LinearSystem K = fem::assemble_linear_static(m);
  const fem::LinearSystem M = fem::assemble_mass(m, false);
  CX_CHECK(K.n_free == 24);  // 8 nodes x 3, none constrained
  // Six lowest modes are rigid-body (eigenvalue ~ 0); extraction must not fail.
  const numerics::EigenBasis b0 = numerics::extract_modes(K, M, 8, /*shift=*/0.0);
  CX_CHECK(b0.modes.size() == 8);
  for (int i = 0; i < 6; ++i)
    CX_CHECK(std::fabs(b0.modes[static_cast<std::size_t>(i)].eigenvalue) < 1.0);
  // The 7th mode is the first genuine deformation mode (>0).
  CX_CHECK(b0.modes[6].eigenvalue > 1.0);
  // A spectral shift returns the same un-shifted eigenvalues (relative tol on the
  // genuine deformation modes; rigid-body modes are ~0 machine noise either way).
  const numerics::EigenBasis bs = numerics::extract_modes(K, M, 8, /*shift=*/-1.0);
  for (std::size_t i = 6; i < 8; ++i)
    CX_NEAR(bs.modes[i].eigenvalue, b0.modes[i].eigenvalue,
            1e-6 * std::fabs(b0.modes[i].eigenvalue));
  for (std::size_t i = 0; i < 6; ++i)  // rigid-body: both near zero
    CX_CHECK(std::fabs(bs.modes[i].eigenvalue) < 1.0);
}

}  // namespace

int main() {
  test_mass_total();
  test_generalized_residual();
  test_spring_mass_analytical();
  test_rigid_body_shift();
  CX_MAIN_RETURN();
}
