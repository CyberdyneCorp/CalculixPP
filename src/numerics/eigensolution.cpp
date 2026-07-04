#include "calculixpp/numerics/eigensolution.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "numpp/core/creation.hpp"
#include "numpp/core/dtype.hpp"
#include "numpp/core/ndarray.hpp"
#include "numpp/linalg/linalg.hpp"
#include "scipp/sparse/sparse.hpp"

// Generalized symmetric eigensolution (spec: eigensolution). Solves K x = λ M x for
// the lowest N modes via SciPP's sparse thick-restart shift-invert Lanczos
// (scipp::sparse::eigsh, SciPP#12): (K - σ M) is factored once (reusing SciPP's sparse
// direct factorization) and applied as the Krylov operator, so only the k lowest modes
// are extracted and the cost/memory stay sparse — scalable past the old dense O(n³)
// path (which relied on NumPP's dense eigh, since fixed to O(n³) in NumPP#138). σ (shift)
// selects the target: 0 for the lowest modes of an SPD pencil, a small negative shift to
// move off a rigid-body / near-zero mode. Eigenvectors come back mass-normalized
// (xᵀ M x = 1) and eigenvalues ascending.
//
// SciPP#15 made eigsh robust for the two FE-pencil failure modes we hit: stiff pencils
// (λ ~ 1e10..1e13 in consistent units, which used to trip an absolute breakdown test)
// now use a breakdown test relative to the operator scale, and tightly clustered spectra
// converge via thick restart — so no manual pencil rescaling is needed. A dense
// M-Cholesky reduction is still kept as a fallback for the one case shift-invert cannot
// handle: a singular (K − σ M) at σ = 0 with rigid-body modes (used when n_free is small
// enough to afford the O(n³) reduction).
namespace cxpp::numerics {
namespace {

// Build a SciPP CSR matrix from a LinearSystem's free/free COO triplets (the same
// full-symmetric COO the sparse solve consumes).
scipp::sparse::CsrMatrix csr_from_system(const fem::LinearSystem& sys) {
  const std::int64_t n = static_cast<std::int64_t>(sys.n_free);
  const std::int64_t nnz = static_cast<std::int64_t>(sys.vals.size());
  scipp::sparse::CooMatrix coo;
  coo.rows = n;
  coo.cols = n;
  coo.data = numpp::zeros({nnz}, numpp::kFloat64);
  coo.row = numpp::zeros({nnz}, numpp::kInt64);
  coo.col = numpp::zeros({nnz}, numpp::kInt64);
  for (std::int64_t i = 0; i < nnz; ++i) {
    const auto k = static_cast<std::size_t>(i);
    coo.data.set_item<double>({i}, sys.vals[k]);
    coo.row.set_item<std::int64_t>({i}, sys.rows[k]);
    coo.col.set_item<std::int64_t>({i}, sys.cols[k]);
  }
  return scipp::sparse::CsrMatrix::from_coo(coo);
}

// Sparse mat-vec y = A x over a LinearSystem's COO (for participation).
std::vector<Real> sparse_matvec(const fem::LinearSystem& A, const std::vector<Real>& x) {
  std::vector<Real> y(static_cast<std::size_t>(A.n_free), 0.0);
  for (std::size_t k = 0; k < A.vals.size(); ++k)
    y[static_cast<std::size_t>(A.rows[k])] +=
        A.vals[k] * x[static_cast<std::size_t>(A.cols[k])];
  return y;
}

// Column j of an (n, k) ndarray into a std::vector<Real>.
std::vector<Real> column(const numpp::ndarray& A, std::int64_t j, std::int64_t n) {
  std::vector<Real> col(static_cast<std::size_t>(n));
  for (std::int64_t i = 0; i < n; ++i)
    col[static_cast<std::size_t>(i)] = A.item<double>({i, j});
  return col;
}

// Expand a free-DOF eigenvector φ to a full nodal mode shape through the constraint
// transform (a mode has zero prescribed displacement; MPC slaves from masters).
std::vector<Vec3> expand_shape(const fem::LinearSystem& sys,
                               const std::vector<Real>& phi, std::size_t n_nodes) {
  const std::vector<Real> zero_prescribed(sys.prescribed.size(), 0.0);
  std::vector<Vec3> shape(n_nodes, Vec3{0, 0, 0});
  for (std::size_t ni = 0; ni < n_nodes; ++ni)
    for (int c = 0; c < kDofsPerNode; ++c) {
      const std::size_t g = ni * kDofsPerNode + static_cast<std::size_t>(c);
      shape[ni][static_cast<std::size_t>(c)] =
          sys.transform.displacement(g, phi, zero_prescribed);
    }
  return shape;
}

// Above this free-DOF count the dense fallback (O(n³) reduction) is too costly, so a
// non-converging sparse eigensolve is reported instead of silently going dense.
constexpr std::int64_t kDenseFallbackMaxDof = 1500;

// Build a dense symmetric (n × n) ndarray from a LinearSystem's free/free COO triplets
// (duplicates summed, both triangles mirrored — the congruence scatter stores some
// entries on a single triangle).
numpp::ndarray dense_from_system(const fem::LinearSystem& sys) {
  const std::int64_t n = static_cast<std::int64_t>(sys.n_free);
  numpp::ndarray A = numpp::zeros({n, n}, numpp::kFloat64);
  for (std::size_t k = 0; k < sys.vals.size(); ++k)
    A.set_item<double>({sys.rows[k], sys.cols[k]},
                       A.item<double>({sys.rows[k], sys.cols[k]}) + sys.vals[k]);
  for (std::int64_t i = 0; i < n; ++i)
    for (std::int64_t j = i + 1; j < n; ++j) {
      const double aij = A.item<double>({i, j}), aji = A.item<double>({j, i});
      const double val = (aij != 0.0 && aji != 0.0) ? aij : aij + aji;
      A.set_item<double>({i, j}, val);
      A.set_item<double>({j, i}, val);
    }
  return A;
}

// Dense generalized symmetric eigensolution K x = λ M x for the lowest num_modes:
//   M = L Lᵀ ; A = L⁻¹ (K − σ M) L⁻ᵀ ; A z = μ z (ascending) ; φ = L⁻ᵀ z.
// M is always SPD, so this handles a singular / rigid-body K that shift-invert Lanczos
// cannot factor. O(n³) — the small-problem fallback and validation oracle.
EigenBasis dense_extract_modes(const fem::LinearSystem& K, const fem::LinearSystem& M,
                               std::size_t num_modes, Real shift) {
  const std::int64_t n = static_cast<std::int64_t>(K.n_free);
  const numpp::ndarray Kd = dense_from_system(K);
  const numpp::ndarray Md = dense_from_system(M);
  numpp::ndarray Ks = Kd;
  if (shift != 0.0)
    for (std::int64_t i = 0; i < n; ++i)
      for (std::int64_t j = 0; j < n; ++j)
        Ks.set_item<double>({i, j},
                            Kd.item<double>({i, j}) - shift * Md.item<double>({i, j}));

  numpp::ndarray L;
  try {
    L = numpp::linalg::cholesky(Md);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("extract_modes: mass matrix is not positive "
                                         "definite (Cholesky failed): ") +
                             e.what());
  }
  // A = L⁻¹ Ks L⁻ᵀ ; back-transform φ = L⁻ᵀ z.
  const numpp::ndarray Y = numpp::linalg::solve(L, Ks);
  const numpp::ndarray A = numpp::linalg::solve(L, Y.transpose()).transpose();
  const numpp::linalg::EighResult eig = numpp::linalg::eigh(A);
  const numpp::ndarray Phi = numpp::linalg::solve(L.transpose(), eig.eigenvectors);

  EigenBasis basis;
  basis.n_free = K.n_free;
  basis.stiffness = K;
  const std::size_t n_nodes = K.dof_eq.size() / static_cast<std::size_t>(kDofsPerNode);
  for (std::size_t m = 0; m < num_modes; ++m) {
    const std::int64_t j = static_cast<std::int64_t>(m);
    std::vector<Real> phi = column(Phi, j, n);
    const std::vector<Real> Mphi = sparse_matvec(M, phi);  // mass-normalize φᵀ M φ = 1
    Real mm = 0.0;
    for (std::size_t i = 0; i < phi.size(); ++i) mm += phi[i] * Mphi[i];
    const Real mnorm = std::sqrt(std::max(mm, 0.0));
    if (mnorm > 0.0)
      for (Real& v : phi) v /= mnorm;

    Mode mode;
    mode.eigenvalue = eig.eigenvalues.item<double>({j}) + shift;  // un-shift
    const Real lam = std::max(mode.eigenvalue, 0.0);
    mode.omega = std::sqrt(lam);
    mode.frequency = mode.omega / (2.0 * M_PI);
    mode.shape = expand_shape(K, phi, n_nodes);
    basis.modes.push_back(std::move(mode));
    basis.mode_free.push_back(std::move(phi));
  }
  return basis;
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
  const int k = static_cast<int>(num_modes);

  // Sparse thick-restart shift-invert Lanczos (SciPP eigsh): lowest k modes of
  // K x = λ M x nearest σ = shift. SciPP#15 made eigsh robust for stiff FE pencils
  // (breakdown test relative to the operator scale) and for tightly clustered spectra
  // (thick restart), so it runs directly on the unscaled pencil — no manual rescaling.
  // maxiter=0 takes SciPP's default restart-cycle budget (the Krylov subspace size is
  // chosen internally); the thick restart converges the k wanted modes across cycles.
  const scipp::sparse::CsrMatrix Kc = csr_from_system(K);
  const scipp::sparse::CsrMatrix Mc = csr_from_system(M);
  const scipp::sparse::EigshResult res =
      scipp::sparse::eigsh(Kc, Mc, k, /*sigma=*/shift, /*tol=*/1e-6, /*maxiter=*/0);
  // (K − σ M) is singular when σ = 0 meets rigid-body modes, which shift-invert cannot
  // factor; fall back to the dense reduction (M is always SPD) for small such problems.
  if (!res.converged) {
    if (n <= kDenseFallbackMaxDof) return dense_extract_modes(K, M, num_modes, shift);
    throw std::runtime_error(
        "extract_modes: eigsh did not converge for the requested " + std::to_string(k) +
        " modes (subspace " + std::to_string(res.iterations) + ") and n_free=" +
        std::to_string(n) + " exceeds the dense-fallback limit; the pencil is likely "
        "singular at this shift — use a shift near the target band");
  }

  EigenBasis basis;
  basis.n_free = K.n_free;
  basis.stiffness = K;
  const std::size_t n_nodes =
      K.dof_eq.size() / static_cast<std::size_t>(kDofsPerNode);

  for (int m = 0; m < k; ++m) {
    Mode mode;
    mode.eigenvalue = res.eigenvalues.item<double>({m});  // eigsh returns λ, ascending
    const Real lam = std::max(mode.eigenvalue, 0.0);      // clamp tiny-negative rigid body
    mode.omega = std::sqrt(lam);
    mode.frequency = mode.omega / (2.0 * M_PI);
    std::vector<Real> phi = column(res.eigenvectors, m, n);  // mass-normalized (φᵀMφ=1)
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

  // Influence vector r: restriction of the rigid-body unit translation on `dir` (at
  // every node) to the free DOFs via the constraint transform (MPC masters accumulate
  // slave terms).
  std::vector<Real> r(static_cast<std::size_t>(n), 0.0);
  const fem::ConstraintTransform& tf = basis.stiffness.transform;
  const std::size_t n_nodes =
      basis.stiffness.dof_eq.size() / static_cast<std::size_t>(kDofsPerNode);
  for (std::size_t ni = 0; ni < n_nodes; ++ni) {
    const std::size_t g = ni * kDofsPerNode + static_cast<std::size_t>(dir);
    if (g >= tf.expansion.size()) continue;
    for (const fem::DofTerm& t : tf.expansion[g].terms)
      r[static_cast<std::size_t>(t.eq)] += t.coeff;
  }
  const std::vector<Real> Mr = sparse_matvec(M, r);  // M r (sparse)

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
