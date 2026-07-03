#pragma once
#include <cstddef>
#include <string>
#include <vector>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/types.hpp"
#include "calculixpp/fem/assembly.hpp"

// Substructure / superelement generation (spec: substructure-generation, tasks 5.1-5.3).
// A *SUBSTRUCTURE GENERATE step condenses the full model onto the retained (master)
// DOFs of *RETAINED NODAL DOFS. Partition the free DOFs into boundary b (retained) and
// interior i (condensed). Static (Guyan) reduction forms the reduced stiffness as the
// Schur complement
//     K_hat = K_bb - K_bi K_ii⁻¹ K_ib,
// and — when a mass matrix is requested — Craig-Bampton reduction builds the constraint
// modes Ψ = -K_ii⁻¹ K_ib and (optionally) fixed-interface normal modes Φ from the
// eigensolution engine, then projects K and M onto the transformation basis
//     T = [ I   0 ]      (retained boundary DOFs + generalized modal DOFs)
//         [ Ψ   Φ ]
// giving the reduced K/M whose leading b×b block reproduces the Guyan stiffness and
// whose modal block is diagonal (ω_k² on K, 1 on M for the mass-normalized modes).
//
// The interior condensation K_ii⁻¹ routes through NumPP (dense factor/solve on the
// ComputeBackend — no GPU needed). This is the DENSE path, O(n_interior³); the same
// scalability note as the eigensolution engine applies (a sparse Schur-complement via
// SciPP is a follow-up). The reduced matrices are exported by write_substructure_matrix
// in the reference lower-triangular *MATRIX TYPE=STIFFNESS/MASS row format.
namespace cxpp::numerics {

// A generated superelement: the reduced operators on the retained DOFs (plus any
// Craig-Bampton generalized modal DOFs appended after the retained block). `n_retained`
// is the number of retained boundary DOFs (the leading block); `n_modes` is the number
// of fixed-interface normal modes appended (0 for pure Guyan). The reduced order is
// `dim = n_retained + n_modes`. `k_reduced` / `m_reduced` are the dim×dim symmetric
// reduced stiffness / mass, stored row-major (m_reduced empty when mass was not
// requested). `retained_node` / `retained_comp` name each retained DOF in output order.
struct Superelement {
  std::size_t n_retained{0};
  std::size_t n_modes{0};
  std::size_t dim() const { return n_retained + n_modes; }
  std::vector<Real> k_reduced;  // dim*dim row-major (symmetric)
  std::vector<Real> m_reduced;  // dim*dim row-major (symmetric); empty if no mass
  std::vector<Index> retained_node;  // node id per retained DOF (size n_retained)
  std::vector<int> retained_comp;    // 1..3 per retained DOF (size n_retained)
  // Fixed-interface modal frequencies ω_k (rad/time) for the appended modal DOFs,
  // size n_modes — the diagonal of the modal stiffness block is ω_k².
  std::vector<Real> modal_omega;
};

// Generate the superelement for a *SUBSTRUCTURE GENERATE model. Requires
// `model.retained_dofs` non-empty (throws std::runtime_error otherwise, per spec —
// "the retained DOF set is empty"). Assembles K (and M when model.substructure_mass),
// resolves each retained (node, comp) to its free equation, Schur-condenses the
// interior block, and — when model.substructure_modes > 0 and mass is present — appends
// `model.substructure_modes` fixed-interface normal modes (Craig-Bampton). A retained
// DOF that is SPC-constrained or is an MPC slave (no free equation) throws.
Superelement generate_substructure(const Model& model);

// Core reduction on already-assembled systems, for testing / reuse. `retained_eq` lists
// the free-equation index of each retained DOF in output order (all distinct, each in
// [0, K.n_free)). `M` may be empty (nullptr) for a stiffness-only Guyan reduction;
// `num_modes` fixed-interface modes are appended only when M is provided. The retained
// node/comp labels are copied from `labels` if non-empty (else left empty).
Superelement reduce_substructure(const fem::LinearSystem& K,
                                 const fem::LinearSystem* M,
                                 const std::vector<Index>& retained_eq,
                                 std::size_t num_modes,
                                 const std::vector<Model::RetainedDof>& labels);

// Write a reduced matrix in the reference CalculiX substructure export format: a
// "*MATRIX TYPE=<type>" header followed by the LOWER triangle row by row (row i holds
// its i+1 entries a[i][0..i]), matching *SUBSTRUCTURE MATRIX OUTPUT (.dat). `type` is
// "STIFFNESS" or "MASS". `mat` is a dim×dim symmetric row-major matrix. (spec:
// substructure-generation — reduced matrix export.)
std::string format_substructure_matrix(const std::vector<Real>& mat, std::size_t dim,
                                       const std::string& type);

// Write the superelement's reduced stiffness (and mass, if present) to `path` in the
// reference format (STIFFNESS block first, then MASS). (spec: substructure-generation —
// *SUBSTRUCTURE MATRIX OUTPUT / write reduced matrices to file.)
void write_substructure_matrix(const std::string& path, const Superelement& se);

}  // namespace cxpp::numerics
