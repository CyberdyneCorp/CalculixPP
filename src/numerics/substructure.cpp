#include "calculixpp/numerics/substructure.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <stdexcept>
#include <unordered_set>

#include "calculixpp/numerics/eigensolution.hpp"
#include "numpp/core/creation.hpp"
#include "numpp/core/dtype.hpp"
#include "numpp/core/ndarray.hpp"
#include "numpp/linalg/linalg.hpp"

// Craig-Bampton / Guyan substructure reduction (spec: substructure-generation). The
// interior condensation and the modal-block projection route through NumPP dense
// factor/solve; the fixed-interface normal modes come from the eigensolution engine
// (extract_modes on the interior-restricted K/M). See the header for the math.
namespace cxpp::numerics {
namespace {

// Dense symmetric matrix (n×n, float64) scattered from a LinearSystem's free/free COO,
// mirrored to both triangles — same construction the eigensolution engine uses so the
// reduced operators match the assembled K/M exactly.
std::vector<Real> dense_symmetric(const fem::LinearSystem& sys) {
  const std::size_t n = static_cast<std::size_t>(sys.n_free);
  std::vector<Real> A(n * n, 0.0);
  for (std::size_t k = 0; k < sys.vals.size(); ++k)
    A[static_cast<std::size_t>(sys.rows[k]) * n +
      static_cast<std::size_t>(sys.cols[k])] += sys.vals[k];
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = i + 1; j < n; ++j) {
      const Real aij = A[i * n + j], aji = A[j * n + i];
      const Real v = (aij != 0.0) ? aij : aji;  // one triangle populated per entry
      A[i * n + j] = A[j * n + i] = v;
    }
  return A;
}

// A dense ndarray (r×c) from a row-major submatrix A[rows_i, cols_j] of the full n×n
// dense matrix `full`. Used to slice K_bb, K_bi, K_ii, K_ib for the Schur complement.
numpp::ndarray sub_ndarray(const std::vector<Real>& full, std::size_t n,
                           const std::vector<std::size_t>& ri,
                           const std::vector<std::size_t>& cj) {
  const std::int64_t r = static_cast<std::int64_t>(ri.size());
  const std::int64_t c = static_cast<std::int64_t>(cj.size());
  numpp::ndarray S = numpp::zeros({r, c}, numpp::kFloat64);
  for (std::int64_t a = 0; a < r; ++a)
    for (std::int64_t b = 0; b < c; ++b)
      S.set_item<double>({a, b}, full[ri[static_cast<std::size_t>(a)] * n +
                                       cj[static_cast<std::size_t>(b)]]);
  return S;
}

// A LinearSystem restricted to a subset of free equations `keep` (the interior DOFs),
// re-numbered 0..|keep|-1, so the eigensolution engine can extract the fixed-interface
// normal modes (K_ii x = λ M_ii x — all interior DOFs "fixed" at the boundary). The
// restricted problem is unconstrained on its own free numbering, so we build a FRESH
// IDENTITY constraint transform whose global DOF space is exactly the |keep| interior
// equations (one node = one DOF, kDofsPerNode padding). This keeps extract_modes'
// shape-expansion self-consistent — we consume the free-DOF eigenvectors (mode_free),
// not the (now-synthetic) nodal shapes. The parent transform must NOT be reused: its
// expansion terms reference the parent free numbering, which is out of range here.
fem::LinearSystem restrict_system(const fem::LinearSystem& sys,
                                  const std::vector<std::size_t>& keep) {
  std::vector<std::int64_t> remap(static_cast<std::size_t>(sys.n_free), -1);
  for (std::size_t k = 0; k < keep.size(); ++k)
    remap[keep[k]] = static_cast<std::int64_t>(k);
  fem::LinearSystem out;
  const Index nk = static_cast<Index>(keep.size());
  out.n_free = nk;
  for (std::size_t k = 0; k < sys.vals.size(); ++k) {
    const std::int64_t r = remap[static_cast<std::size_t>(sys.rows[k])];
    const std::int64_t c = remap[static_cast<std::size_t>(sys.cols[k])];
    if (r < 0 || c < 0) continue;  // drop rows/cols touching the boundary block
    out.rows.push_back(static_cast<Index>(r));
    out.cols.push_back(static_cast<Index>(c));
    out.vals.push_back(sys.vals[k]);
  }
  // Fresh identity transform over `nk` global DOFs (each maps to its own free equation).
  // Padded to a whole number of pseudo-nodes so dof_eq.size() % kDofsPerNode == 0 (the
  // eigensolution shape-expansion divides by kDofsPerNode). The padding DOFs are
  // constrained (dof_eq = -1), contributing nothing.
  const std::size_t ndof = (static_cast<std::size_t>(nk) +
                            static_cast<std::size_t>(kDofsPerNode) - 1) /
                           static_cast<std::size_t>(kDofsPerNode) *
                           static_cast<std::size_t>(kDofsPerNode);
  out.dof_eq.assign(ndof, -1);
  out.prescribed.assign(ndof, 0.0);
  fem::ConstraintTransform& tf = out.transform;
  tf.n_free = nk;
  tf.dof_eq.assign(ndof, -1);
  tf.expansion.assign(ndof, fem::DofExpansion{});
  tf.is_slave.assign(ndof, false);
  tf.spc_terms.assign(ndof, {});
  for (Index e = 0; e < nk; ++e) {
    const std::size_t g = static_cast<std::size_t>(e);
    out.dof_eq[g] = e;
    tf.dof_eq[g] = e;
    tf.expansion[g].terms.push_back(fem::DofTerm{e, 1.0});
  }
  return out;
}

// Symmetrize a dim×dim row-major matrix in place (average off-diagonals) so tiny
// asymmetry from the dense solve does not leak into the exported reduced matrix.
void symmetrize(std::vector<Real>& A, std::size_t dim) {
  for (std::size_t i = 0; i < dim; ++i)
    for (std::size_t j = i + 1; j < dim; ++j) {
      const Real s = 0.5 * (A[i * dim + j] + A[j * dim + i]);
      A[i * dim + j] = A[j * dim + i] = s;
    }
}

// Elementwise A - B for two same-shape (r×c) ndarrays.
numpp::ndarray sub_elem(const numpp::ndarray& A, const numpp::ndarray& B,
                        std::int64_t r, std::int64_t c) {
  numpp::ndarray R = numpp::zeros({r, c}, numpp::kFloat64);
  for (std::int64_t i = 0; i < r; ++i)
    for (std::int64_t j = 0; j < c; ++j)
      R.set_item<double>({i, j}, A.item<double>({i, j}) - B.item<double>({i, j}));
  return R;
}

// Copy the leading nb×nb block of an ndarray into the top-left of a dim×dim row-major
// buffer (the retained/retained block of the reduced K or M).
void scatter_block(std::vector<Real>& dst, std::size_t dim, const numpp::ndarray& A,
                   std::size_t nb) {
  for (std::size_t i = 0; i < nb; ++i)
    for (std::size_t j = 0; j < nb; ++j)
      dst[i * dim + j] = A.item<double>({static_cast<std::int64_t>(i),
                                         static_cast<std::int64_t>(j)});
}

// The interior partition of a Guyan/Craig-Bampton reduction: the reduced-stiffness Schur
// complement K̂ = K_bb − K_bi K_ii⁻¹ K_ib and the constraint modes Ψ = −K_ii⁻¹ K_ib
// (n_i × nb), from the dense K sliced on boundary set `b` and interior set `i`.
struct InteriorReduction {
  numpp::ndarray khat;  // nb × nb reduced stiffness
  numpp::ndarray psi;   // n_i × nb constraint modes
};
InteriorReduction guyan_schur(const std::vector<Real>& Kd, std::size_t n,
                              const std::vector<std::size_t>& b,
                              const std::vector<std::size_t>& interior) {
  const std::int64_t nb = static_cast<std::int64_t>(b.size());
  const std::int64_t ni = static_cast<std::int64_t>(interior.size());
  const numpp::ndarray Kbb = sub_ndarray(Kd, n, b, b);
  const numpp::ndarray Kbi = sub_ndarray(Kd, n, b, interior);
  const numpp::ndarray Kib = sub_ndarray(Kd, n, interior, b);
  const numpp::ndarray Kii = sub_ndarray(Kd, n, interior, interior);
  const numpp::ndarray KiiInvKib = numpp::linalg::solve(Kii, Kib);  // K_ii⁻¹ K_ib
  InteriorReduction r;
  r.khat = sub_elem(Kbb, numpp::dot(Kbi, KiiInvKib), nb, nb);        // K_bb − K_bi ...
  r.psi = numpp::zeros({ni, nb}, numpp::kFloat64);
  for (std::int64_t i = 0; i < ni; ++i)
    for (std::int64_t j = 0; j < nb; ++j)
      r.psi.set_item<double>({i, j}, -KiiInvKib.item<double>({i, j}));  // Ψ = −K_ii⁻¹K_ib
  return r;
}

// Project the mass matrix onto the Craig-Bampton basis T = [[I,0],[Ψ,Φ]] and write the
// dim×dim reduced mass (row-major) into `Mr`. `Md` is the dense full M; `b`/`interior`
// the partition; `Psi` the constraint modes (empty ndarray when no interior); `phi_i` the
// fixed-interface mode vectors over the interior DOFs (mass-normalized). Blocks:
//   M̂_bb = M_bb + M_bi Ψ + Ψᵀ M_ib + Ψᵀ M_ii Ψ, M_bm = (M_bi + Ψᵀ M_ii) Φ, M_mm = I.
void project_reduced_mass(std::vector<Real>& Mr, std::size_t dim,
                          const std::vector<Real>& Md, std::size_t n,
                          const std::vector<std::size_t>& b,
                          const std::vector<std::size_t>& interior,
                          const numpp::ndarray& Psi,
                          const std::vector<std::vector<Real>>& phi_i) {
  const std::size_t nb = b.size();
  const std::int64_t NB = static_cast<std::int64_t>(nb);
  const std::int64_t NI = static_cast<std::int64_t>(interior.size());
  const numpp::ndarray Mbb = sub_ndarray(Md, n, b, b);
  numpp::ndarray Mhat_bb = Mbb;
  numpp::ndarray Mib_col;  // (M_bi + Ψᵀ M_ii)  nb × n_i (modal coupling row)
  if (!interior.empty()) {
    const numpp::ndarray Mbi = sub_ndarray(Md, n, b, interior);
    const numpp::ndarray Mib = sub_ndarray(Md, n, interior, b);
    const numpp::ndarray Mii = sub_ndarray(Md, n, interior, interior);
    const numpp::ndarray PsiT = Psi.transpose();
    // M̂_bb = M_bb + M_bi Ψ + Ψᵀ M_ib + Ψᵀ M_ii Ψ (four same-shape nb×nb terms).
    Mhat_bb = numpp::zeros({NB, NB}, numpp::kFloat64);
    const numpp::ndarray t1 = numpp::dot(Mbi, Psi);
    const numpp::ndarray t2 = numpp::dot(PsiT, Mib);
    const numpp::ndarray t3 = numpp::dot(PsiT, numpp::dot(Mii, Psi));
    for (std::int64_t i = 0; i < NB; ++i)
      for (std::int64_t j = 0; j < NB; ++j)
        Mhat_bb.set_item<double>({i, j}, Mbb.item<double>({i, j}) +
                                             t1.item<double>({i, j}) +
                                             t2.item<double>({i, j}) +
                                             t3.item<double>({i, j}));
    // M_bi + Ψᵀ M_ii  (nb × n_i).
    const numpp::ndarray PsiTMii = numpp::dot(PsiT, Mii);
    Mib_col = numpp::zeros({NB, NI}, numpp::kFloat64);
    for (std::int64_t i = 0; i < NB; ++i)
      for (std::int64_t j = 0; j < NI; ++j)
        Mib_col.set_item<double>({i, j},
                                 Mbi.item<double>({i, j}) + PsiTMii.item<double>({i, j}));
  }
  scatter_block(Mr, dim, Mhat_bb, nb);
  // Retained/modal coupling M_bm = (M_bi + Ψᵀ M_ii) φ_k, modal/modal = 1.
  for (std::size_t k = 0; k < phi_i.size(); ++k) {
    const std::size_t d = nb + k;
    const std::vector<Real>& phik = phi_i[k];
    for (std::size_t i = 0; i < nb; ++i) {
      Real c = 0.0;
      for (std::int64_t j = 0; j < NI; ++j)
        c += Mib_col.item<double>({static_cast<std::int64_t>(i), j}) *
             phik[static_cast<std::size_t>(j)];
      Mr[i * dim + d] = c;
      Mr[d * dim + i] = c;
    }
    Mr[d * dim + d] = 1.0;  // φ_kᵀ M_ii φ_k = 1 (mass-normalized)
  }
  symmetrize(Mr, dim);
}

// Partition the free equations into the boundary (retained) set `b` — the deduplicated,
// in-range retained equations in declaration order — and the interior complement.
// Throws on an out-of-range or duplicate retained equation.
void partition_dofs(const std::vector<Index>& retained_eq, std::size_t n,
                    std::vector<std::size_t>& b, std::vector<std::size_t>& interior) {
  std::unordered_set<std::size_t> in_b;
  b.reserve(retained_eq.size());
  for (const Index e : retained_eq) {
    const std::size_t eq = static_cast<std::size_t>(e);
    if (eq >= n) throw std::runtime_error("reduce_substructure: retained eq out of range");
    if (!in_b.insert(eq).second)
      throw std::runtime_error("reduce_substructure: duplicate retained DOF");
    b.push_back(eq);
  }
  interior.reserve(n - b.size());
  for (std::size_t e = 0; e < n; ++e)
    if (in_b.find(e) == in_b.end()) interior.push_back(e);
}

// Extract the fixed-interface normal modes (Craig-Bampton): the lowest `num_modes` modes
// of the interior-restricted K_ii x = λ M_ii x. Fills `omega` / `phi_i` (mode vectors on
// the interior DOFs) and returns the count actually extracted.
std::size_t fixed_interface_modes(const fem::LinearSystem& K, const fem::LinearSystem& M,
                                  const std::vector<std::size_t>& interior,
                                  std::size_t num_modes, std::vector<Real>& omega,
                                  std::vector<std::vector<Real>>& phi_i) {
  const fem::LinearSystem Ki = restrict_system(K, interior);
  const fem::LinearSystem Mi = restrict_system(M, interior);
  const EigenBasis fib =
      extract_modes(Ki, Mi, std::min(num_modes, interior.size()));
  for (const Mode& mode : fib.modes) omega.push_back(mode.omega);
  phi_i = fib.mode_free;
  return fib.modes.size();
}

}  // namespace

Superelement reduce_substructure(const fem::LinearSystem& K,
                                 const fem::LinearSystem* M,
                                 const std::vector<Index>& retained_eq,
                                 std::size_t num_modes,
                                 const std::vector<Model::RetainedDof>& labels) {
  const std::size_t n = static_cast<std::size_t>(K.n_free);
  const std::size_t nb = retained_eq.size();
  if (nb == 0) throw std::runtime_error("reduce_substructure: empty retained DOF set");

  std::vector<std::size_t> b, interior;
  partition_dofs(retained_eq, n, b, interior);
  const bool have_interior = !interior.empty();

  const std::vector<Real> Kd = dense_symmetric(K);
  // Guyan Schur complement K̂ + constraint modes Ψ (identity K̂ = K_bb, empty Ψ when
  // there is no interior block).
  numpp::ndarray Khat = sub_ndarray(Kd, n, b, b);
  numpp::ndarray Psi;
  if (have_interior) {
    const InteriorReduction ir = guyan_schur(Kd, n, b, interior);
    Khat = ir.khat;
    Psi = ir.psi;
  }

  Superelement se;
  se.n_retained = nb;
  for (std::size_t k = 0; k < nb && k < labels.size(); ++k) {
    se.retained_node.push_back(labels[k].node_id);
    se.retained_comp.push_back(labels[k].comp);
  }

  // Fixed-interface normal modes (Craig-Bampton) — only with a mass matrix, an interior
  // block, and a nonzero mode request.
  std::vector<std::vector<Real>> phi_i;  // per mode, size n_interior (mass-normalized)
  if (M != nullptr && have_interior && num_modes > 0)
    se.n_modes = fixed_interface_modes(K, *M, interior, num_modes, se.modal_omega, phi_i);

  const std::size_t dim = se.dim();
  // Reduced stiffness: leading nb×nb Guyan block + diagonal ω_k² modal block (zero
  // coupling — Craig-Bampton block-diagonalizes K by construction).
  std::vector<Real> Kr(dim * dim, 0.0);
  scatter_block(Kr, dim, Khat, nb);
  for (std::size_t k = 0; k < se.n_modes; ++k)
    Kr[(nb + k) * dim + (nb + k)] = se.modal_omega[k] * se.modal_omega[k];  // ω_k²
  symmetrize(Kr, dim);
  se.k_reduced = std::move(Kr);

  // Reduced mass (Craig-Bampton) — project M onto T = [[I,0],[Ψ,Φ]].
  if (M != nullptr) {
    std::vector<Real> Mr(dim * dim, 0.0);
    project_reduced_mass(Mr, dim, dense_symmetric(*M), n, b, interior, Psi, phi_i);
    se.m_reduced = std::move(Mr);
  }

  return se;
}

Superelement generate_substructure(const Model& model) {
  if (model.retained_dofs.empty())
    throw std::runtime_error(
        "*SUBSTRUCTURE GENERATE: the retained DOF set is empty (declare "
        "*RETAINED NODAL DOFS)");

  const fem::LinearSystem K = fem::assemble_linear_static(model);
  fem::LinearSystem M;
  const bool want_mass = model.substructure_mass || model.substructure_modes > 0;
  if (want_mass) M = fem::assemble_mass(model, /*lumped=*/false);

  // Resolve each retained (node, comp) to its free equation via the K DOF map.
  std::vector<Index> retained_eq;
  retained_eq.reserve(model.retained_dofs.size());
  for (const Model::RetainedDof& rd : model.retained_dofs) {
    const std::size_t ni = static_cast<std::size_t>(model.mesh.node_index(rd.node_id));
    const std::size_t g =
        ni * static_cast<std::size_t>(kDofsPerNode) + static_cast<std::size_t>(rd.comp - 1);
    if (g >= K.dof_eq.size() || K.dof_eq[g] < 0)
      throw std::runtime_error(
          "*RETAINED NODAL DOFS: retained DOF is constrained (SPC) or an MPC slave — "
          "it must remain a free (active) DOF of the superelement (node " +
          std::to_string(rd.node_id) + ", dof " + std::to_string(rd.comp) + ")");
    retained_eq.push_back(K.dof_eq[g]);
  }

  return reduce_substructure(K, want_mass ? &M : nullptr, retained_eq,
                             static_cast<std::size_t>(std::max(0, model.substructure_modes)),
                             model.retained_dofs);
}

std::string format_substructure_matrix(const std::vector<Real>& mat, std::size_t dim,
                                       const std::string& type) {
  std::string out = "*MATRIX TYPE=" + type + "\n";
  char buf[32];
  for (std::size_t i = 0; i < dim; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      std::snprintf(buf, sizeof(buf), " % .13E", mat[i * dim + j]);
      out += buf;
    }
    out += "\n";
  }
  return out;
}

void write_substructure_matrix(const std::string& path, const Superelement& se) {
  std::ofstream f(path);
  if (!f) throw std::runtime_error("write_substructure_matrix: cannot open " + path);
  const std::size_t dim = se.dim();
  f << format_substructure_matrix(se.k_reduced, dim, "STIFFNESS");
  if (!se.m_reduced.empty())
    f << format_substructure_matrix(se.m_reduced, dim, "MASS");
}

}  // namespace cxpp::numerics
