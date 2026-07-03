#include "calculixpp/numerics/eigensolution.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>
#include <stdexcept>

#include "numpp/core/creation.hpp"
#include "numpp/core/dtype.hpp"
#include "numpp/core/ndarray.hpp"
#include "numpp/linalg/linalg.hpp"

// Dense generalized symmetric eigensolution (spec: eigensolution). We solve
// K x = λ M x for the lowest N modes via the standard reduction:
//   M = L Lᵀ  (Cholesky, numpp::linalg::cholesky)
//   A = L⁻¹ K L⁻ᵀ  (symmetric standard problem, via triangular solves)
//   A z = μ z  (numpp::linalg::eigh — ascending eigenvalues)
//   x = L⁻ᵀ z  (back-transform)  then mass-normalize (xᵀ M x = 1).
// A spectral shift σ solves (K - σ M) so rigid-body / near-zero (and, for buckling,
// non-positive) eigenvalues are handled; eigenvalues are returned un-shifted.
//
// SCALABILITY [~]: this dense path is O(n_free³) — correct and the validation
// oracle, but not for large meshes. The scalable target is a shift-invert Lanczos
// driving SciPP's sparse factorization ((K - σ M)⁻¹ via scipp::sparse::spsolve) as
// the Krylov operator; NumPP/SciPP expose no sparse GENERALIZED eigensolver yet, so
// that path is a SciPP follow-up (SciPP#12; the dense eigh/solve it would replace
// are O(n^4) — NumPP#138). The dense path below is what *FREQUENCY
// uses today and is validated against a stock CalculiX .dat.ref.
namespace cxpp::numerics {
namespace {

// Build a dense symmetric ndarray (n x n, float64) from the free/free COO triplets
// of a LinearSystem (duplicates summed, both triangles filled). The assembler emits
// each entry once (upper OR lower via the transform), so we symmetrize by summing
// (r,c) and (c,r) into both slots only when they are distinct keys — the COO already
// carries the exact scattered values, so a plain scatter into a dense buffer with an
// explicit lower/upper mirror reproduces the full symmetric matrix.
numpp::ndarray dense_from_system(const fem::LinearSystem& sys) {
  const std::int64_t n = static_cast<std::int64_t>(sys.n_free);
  numpp::ndarray A = numpp::zeros({n, n}, numpp::kFloat64);
  for (std::size_t k = 0; k < sys.vals.size(); ++k) {
    const std::int64_t r = static_cast<std::int64_t>(sys.rows[k]);
    const std::int64_t c = static_cast<std::int64_t>(sys.cols[k]);
    const Real v = sys.vals[k];
    A.set_item<double>({r, c}, A.item<double>({r, c}) + v);
  }
  // Mirror to guarantee exact symmetry (the congruence scatter can leave the matrix
  // stored on one triangle for some entries). A_sym = (A + Aᵀ) - diag(A).
  for (std::int64_t i = 0; i < n; ++i)
    for (std::int64_t j = i + 1; j < n; ++j) {
      const double s = A.item<double>({i, j}) + A.item<double>({j, i});
      // If only one side was populated the other is 0, so the sum is the true value;
      // if both were populated (shouldn't happen for a single-scatter) we'd double —
      // guard by taking the larger-magnitude convention: prefer the nonzero side.
      const double aij = A.item<double>({i, j}), aji = A.item<double>({j, i});
      const double val = (aij != 0.0 && aji != 0.0) ? aij : s;
      A.set_item<double>({i, j}, val);
      A.set_item<double>({j, i}, val);
    }
  return A;
}

// K_shift = K - shift * M  (dense).
numpp::ndarray shifted_stiffness(const numpp::ndarray& K, const numpp::ndarray& M,
                                 Real shift, std::int64_t n) {
  if (shift == 0.0) return K;
  numpp::ndarray Ks = numpp::zeros({n, n}, numpp::kFloat64);
  for (std::int64_t i = 0; i < n; ++i)
    for (std::int64_t j = 0; j < n; ++j)
      Ks.set_item<double>({i, j},
                          K.item<double>({i, j}) - shift * M.item<double>({i, j}));
  return Ks;
}

// Read column j of a matrix (n rows) into a std::vector<Real>.
std::vector<Real> column(const numpp::ndarray& A, std::int64_t j, std::int64_t n) {
  std::vector<Real> col(static_cast<std::size_t>(n));
  for (std::int64_t i = 0; i < n; ++i)
    col[static_cast<std::size_t>(i)] = A.item<double>({i, j});
  return col;
}

// y = M * x  (dense mat-vec on the free DOFs), for mass normalization / participation.
std::vector<Real> mass_matvec(const numpp::ndarray& M, const std::vector<Real>& x,
                              std::int64_t n) {
  std::vector<Real> y(static_cast<std::size_t>(n), 0.0);
  for (std::int64_t i = 0; i < n; ++i) {
    double s = 0.0;
    for (std::int64_t j = 0; j < n; ++j)
      s += M.item<double>({i, j}) * x[static_cast<std::size_t>(j)];
    y[static_cast<std::size_t>(i)] = s;
  }
  return y;
}

Real dot(const std::vector<Real>& a, const std::vector<Real>& b) {
  Real s = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
  return s;
}

// Expand a free-DOF eigenvector φ to a full nodal mode shape through the constraint
// transform (SPC value 0 in a modal problem — a mode has no prescribed displacement;
// MPC slaves reconstructed from masters). Mirrors StaticFields expansion.
std::vector<Vec3> expand_shape(const fem::LinearSystem& sys,
                               const std::vector<Real>& phi, std::size_t n_nodes) {
  std::vector<Real> zero_prescribed(sys.prescribed.size(), 0.0);
  std::vector<Vec3> shape(n_nodes, Vec3{0, 0, 0});
  for (std::size_t ni = 0; ni < n_nodes; ++ni)
    for (int c = 0; c < kDofsPerNode; ++c) {
      const std::size_t g = ni * kDofsPerNode + static_cast<std::size_t>(c);
      shape[ni][static_cast<std::size_t>(c)] =
          sys.transform.displacement(g, phi, zero_prescribed);
    }
  return shape;
}

}  // namespace

EigenBasis extract_modes(const fem::LinearSystem& K, const fem::LinearSystem& M,
                         std::size_t num_modes, Real shift) {
  if (K.n_free != M.n_free)
    throw std::runtime_error(
        "extract_modes: K and M have different free-DOF counts (assemble both from "
        "the same model)");
  const std::int64_t n = static_cast<std::int64_t>(K.n_free);
  if (n == 0) throw std::runtime_error("extract_modes: no free DOFs");
  if (num_modes == 0 || static_cast<std::int64_t>(num_modes) > n)
    num_modes = static_cast<std::size_t>(n);

  const numpp::ndarray Kd = dense_from_system(K);
  const numpp::ndarray Md = dense_from_system(M);
  const numpp::ndarray Ks = shifted_stiffness(Kd, Md, shift, n);

  // M = L Lᵀ ; throws (via numpp) if M is not positive definite.
  numpp::ndarray L;
  try {
    L = numpp::linalg::cholesky(Md);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("extract_modes: mass matrix is not "
                                         "positive definite (Cholesky failed): ") +
                             e.what());
  }
  const numpp::ndarray Lt = L.transpose();

  // A = L⁻¹ Ks L⁻ᵀ. Y = L⁻¹ Ks = solve(L, Ks); A = Y L⁻ᵀ = solve(L, Yᵀ)ᵀ (symmetric).
  const numpp::ndarray Y = numpp::linalg::solve(L, Ks);
  const numpp::ndarray A =
      numpp::linalg::solve(L, Y.transpose()).transpose();

  // Standard symmetric eigenproblem, eigenvalues ascending.
  const numpp::linalg::EighResult eig = numpp::linalg::eigh(A);

  // Back-transform all eigenvectors at once: Z columns z_k ; φ = L⁻ᵀ z = solve(Lᵀ, Z).
  const numpp::ndarray Phi = numpp::linalg::solve(Lt, eig.eigenvectors);

  EigenBasis basis;
  basis.n_free = K.n_free;
  basis.stiffness = K;
  const std::size_t n_nodes =
      K.dof_eq.size() / static_cast<std::size_t>(kDofsPerNode);

  for (std::size_t m = 0; m < num_modes; ++m) {
    const std::int64_t j = static_cast<std::int64_t>(m);
    std::vector<Real> phi = column(Phi, j, n);
    // Mass-normalize: φ <- φ / sqrt(φᵀ M φ). (Un-shifted M — normalization is scale.)
    const std::vector<Real> Mphi = mass_matvec(Md, phi, n);
    const Real mnorm = std::sqrt(std::max(dot(phi, Mphi), 0.0));
    if (mnorm > 0.0)
      for (Real& v : phi) v /= mnorm;

    Mode mode;
    mode.eigenvalue = eig.eigenvalues.item<double>({j}) + shift;  // un-shift
    const Real lam = std::max(mode.eigenvalue, 0.0);  // clamp tiny-negative rigid body
    mode.omega = std::sqrt(lam);
    mode.frequency = mode.omega / (2.0 * M_PI);
    mode.shape = expand_shape(K, phi, n_nodes);

    basis.modes.push_back(std::move(mode));
    basis.mode_free.push_back(std::move(phi));
  }
  return basis;
}

EigenBasis extract_modes(const Model& model, std::size_t num_modes, Real shift) {
  const fem::LinearSystem K = fem::assemble_linear_static(model);
  const fem::LinearSystem M = fem::assemble_mass(model, /*lumped=*/false);
  return extract_modes(K, M, num_modes, shift);
}

Participation participation(const EigenBasis& basis, const fem::LinearSystem& M,
                            int dir) {
  Participation out;
  out.direction = dir;
  const std::int64_t n = static_cast<std::int64_t>(basis.n_free);
  const numpp::ndarray Md = dense_from_system(M);
  const std::size_t n_nodes =
      basis.stiffness.dof_eq.size() / static_cast<std::size_t>(kDofsPerNode);

  // Influence vector r = free-DOF restriction of the rigid-body direction d (unit
  // translation on `dir` at every node). Build d on the full DOF space, then project
  // to the free DOFs via the transform's Tᵀ action: r_eq = Σ_g T[g][eq] d[g]. We
  // realize this by expanding each free equation's unit basis and dotting with d —
  // equivalently, r = restriction of d to the free numbering. For SPC/MPC-free DOFs
  // that is just d[g] at the DOF's own equation; MPC masters accumulate slave terms.
  std::vector<Real> r(static_cast<std::size_t>(n), 0.0);
  const fem::ConstraintTransform& tf = basis.stiffness.transform;
  for (std::size_t ni = 0; ni < n_nodes; ++ni) {
    const std::size_t g = ni * kDofsPerNode + static_cast<std::size_t>(dir);
    if (g >= tf.expansion.size()) continue;
    for (const fem::DofTerm& t : tf.expansion[g].terms)
      r[static_cast<std::size_t>(t.eq)] += t.coeff;  // d[g] = 1 on this direction
  }
  const std::vector<Real> Mr = [&] {
    std::vector<Real> y(static_cast<std::size_t>(n), 0.0);
    for (std::int64_t i = 0; i < n; ++i) {
      double s = 0.0;
      for (std::int64_t j = 0; j < n; ++j)
        s += Md.item<double>({i, j}) * r[static_cast<std::size_t>(j)];
      y[static_cast<std::size_t>(i)] = s;
    }
    return y;
  }();

  out.factor.reserve(basis.modes.size());
  out.effective_mass.reserve(basis.modes.size());
  for (const std::vector<Real>& phi : basis.mode_free) {
    Real g = 0.0;
    for (std::size_t i = 0; i < phi.size(); ++i) g += phi[i] * Mr[i];  // Γ = φᵀ M r
    out.factor.push_back(g);
    out.effective_mass.push_back(g * g);  // mass-normalized -> m_eff = Γ²
    out.total_effective_mass += g * g;
  }
  return out;
}

}  // namespace cxpp::numerics
