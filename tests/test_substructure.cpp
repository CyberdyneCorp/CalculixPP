// Substructure / superelement generation (Phase 4, tasks 5.1-5.3).
//   * Guyan (static) reduction Schur complement against a hand-solvable 3-DOF chain;
//   * empty-retained-DOF-set error (spec scenario);
//   * exported *MATRIX TYPE=STIFFNESS lower-triangular row format;
//   * Craig-Bampton reduced eigenfrequencies approximate the FULL-model *FREQUENCY
//     (the key physics-fidelity check: retaining interface DOFs + fixed-interface
//     normal modes reproduces the low modes of the un-reduced model);
//   * the stock-CalculiX substructure.inp reduced stiffness matches its .dat.ref
//     (validated when the CalculiX test tree is present; skipped gracefully otherwise —
//     the Python regression suite also pins it).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "calculixpp/fem/assembly.hpp"
#include "calculixpp/io/inp_parser.hpp"
#include "calculixpp/numerics/eigensolution.hpp"
#include "calculixpp/numerics/substructure.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

// Build a LinearSystem directly from a dense symmetric matrix (n×n row-major), with an
// identity DOF map (every DOF free) — a minimal harness for reduce_substructure.
fem::LinearSystem system_from_dense(const std::vector<Real>& A, std::size_t n) {
  fem::LinearSystem sys;
  sys.n_free = static_cast<Index>(n);
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j)
      if (A[i * n + j] != 0.0) {
        sys.rows.push_back(static_cast<Index>(i));
        sys.cols.push_back(static_cast<Index>(j));
        sys.vals.push_back(A[i * n + j]);
      }
  return sys;
}

// Guyan reduction of a 3-DOF grounded spring chain, retaining DOFs {0, 2} and
// condensing the middle DOF 1. Springs k between (ground,0),(0,1),(1,2) give
//   K = [[2k,-k,0],[-k,2k,-k],[0,-k,2k]].
// The Schur complement onto {0,2} is K_hat = K_bb - K_bi K_ii⁻¹ K_ib with
//   K_bb=[[2k,0],[0,2k]], K_bi=[[-k],[-k]], K_ii=[2k], K_ib=[[-k],[-k]].
//   K_hat = [[2k,0],[0,2k]] - (1/2k)[[k²,k²],[k²,k²]]
//         = [[2k - k/2, -k/2],[-k/2, 2k - k/2]] = [[1.5k,-0.5k],[-0.5k,1.5k]].
void test_guyan_3dof_analytical() {
  const Real k = 100.0;
  const std::vector<Real> K = {2 * k, -k, 0, -k, 2 * k, -k, 0, -k, 2 * k};
  const fem::LinearSystem sys = system_from_dense(K, 3);
  const numerics::Superelement se =
      numerics::reduce_substructure(sys, nullptr, {0, 2}, 0, {});
  CX_CHECK(se.dim() == 2);
  CX_CHECK(se.n_modes == 0);
  CX_NEAR(se.k_reduced[0], 1.5 * k, 1e-9);
  CX_NEAR(se.k_reduced[1], -0.5 * k, 1e-9);
  CX_NEAR(se.k_reduced[2], -0.5 * k, 1e-9);
  CX_NEAR(se.k_reduced[3], 1.5 * k, 1e-9);
}

// A retained set that covers ALL DOFs (no interior) reduces to the identity: K_hat = K.
void test_no_interior_is_identity() {
  const Real k = 3.0;
  const std::vector<Real> K = {2 * k, -k, -k, 2 * k};
  const fem::LinearSystem sys = system_from_dense(K, 2);
  const numerics::Superelement se =
      numerics::reduce_substructure(sys, nullptr, {0, 1}, 0, {});
  CX_NEAR(se.k_reduced[0], 2 * k, 1e-12);
  CX_NEAR(se.k_reduced[1], -k, 1e-12);
  CX_NEAR(se.k_reduced[3], 2 * k, 1e-12);
}

// Empty retained-DOF set must throw (spec scenario: no retained DOFs declared).
void test_empty_retained_throws() {
  Model m;  // no retained_dofs, procedure Substructure
  m.mesh.add_node(1, {0, 0, 0});
  m.procedure = Procedure::Substructure;
  bool threw = false;
  try {
    numerics::generate_substructure(m);
  } catch (const std::exception&) {
    threw = true;
  }
  CX_CHECK(threw);
}

// Exported reduced-stiffness format: lower triangle, row i has i+1 entries.
void test_export_format_lower_triangle() {
  numerics::Superelement se;
  se.n_retained = 2;
  se.n_modes = 0;
  se.k_reduced = {1.0, 2.0, 2.0, 4.0};  // symmetric 2×2
  const std::string txt =
      numerics::format_substructure_matrix(se.k_reduced, 2, "STIFFNESS");
  CX_CHECK(txt.find("*MATRIX TYPE=STIFFNESS") != std::string::npos);
  // Parse the numeric rows: row 0 -> 1 value, row 1 -> 2 values.
  std::istringstream is(txt);
  std::string header;
  std::getline(is, header);
  std::vector<std::vector<Real>> rows;
  std::string line;
  while (std::getline(is, line)) {
    if (line.empty()) continue;
    std::istringstream ls(line);
    std::vector<Real> vals;
    Real v;
    while (ls >> v) vals.push_back(v);
    if (!vals.empty()) rows.push_back(vals);
  }
  CX_CHECK(rows.size() == 2);
  CX_CHECK(rows[0].size() == 1 && rows[1].size() == 2);
  CX_NEAR(rows[0][0], 1.0, 1e-12);
  CX_NEAR(rows[1][0], 2.0, 1e-12);
  CX_NEAR(rows[1][1], 4.0, 1e-12);
}

// Build a small cantilever C3D8 block (2 stacked hexes), extract the FULL model's
// natural frequencies, then Craig-Bampton reduce onto the free top-face DOFs plus a
// handful of fixed-interface modes and confirm the reduced eigenfrequencies approximate
// the full-model low frequencies. This is the physics-fidelity check for tasks 5.1-5.2.
void test_craig_bampton_approximates_full() {
  Model m;
  // 12 nodes: 3 layers (z=0 clamped, z=1, z=2) of a unit-square column.
  const Real dz = 1.0;
  int id = 1;
  auto add_layer = [&](Real z) {
    m.mesh.add_node(id++, {0, 0, z});
    m.mesh.add_node(id++, {1, 0, z});
    m.mesh.add_node(id++, {1, 1, z});
    m.mesh.add_node(id++, {0, 1, z});
  };
  add_layer(0.0);
  add_layer(dz);
  add_layer(2 * dz);
  // Two stacked C3D8 hexes.
  m.mesh.add_element(1, ElementType::C3D8, {1, 2, 3, 4, 5, 6, 7, 8});
  m.mesh.add_element(2, ElementType::C3D8, {5, 6, 7, 8, 9, 10, 11, 12});
  m.mesh.add_elset("EALL", {1, 2});
  Material mat;
  mat.name = "EL";
  mat.elastic = ElasticIso{210000.0, 0.3};
  mat.density = 7.8e-9;
  m.materials["EL"] = mat;
  m.sections.push_back(SolidSection{"EALL", "EL"});
  // Clamp the z=0 face (nodes 1..4).
  for (Index n : {Index{1}, Index{2}, Index{3}, Index{4}})
    for (int c = 1; c <= 3; ++c) m.spcs.push_back(Spc{n, c, 0.0});
  m.procedure = Procedure::Frequency;

  const fem::LinearSystem K = fem::assemble_linear_static(m);
  const fem::LinearSystem M = fem::assemble_mass(m, /*lumped=*/false);
  const numerics::EigenBasis full = numerics::extract_modes(K, M, 6);

  // Retain the TOP-face translational DOFs (nodes 9..12, dofs 1..3) = 12 boundary DOFs,
  // plus 6 fixed-interface normal modes -> Craig-Bampton reduced order 18 (full free =
  // 8 nodes*3 = 24). Resolve retained DOFs to free equations via the K DOF map.
  std::vector<Index> retained_eq;
  std::vector<Model::RetainedDof> labels;
  for (Index nid : {Index{9}, Index{10}, Index{11}, Index{12}})
    for (int c = 1; c <= 3; ++c) {
      const std::size_t g =
          static_cast<std::size_t>(m.mesh.node_index(nid)) * 3 +
          static_cast<std::size_t>(c - 1);
      retained_eq.push_back(K.dof_eq[g]);
      labels.push_back(Model::RetainedDof{nid, c});
    }
  const numerics::Superelement se =
      numerics::reduce_substructure(K, &M, retained_eq, 6, labels);
  CX_CHECK(se.n_retained == 12);
  CX_CHECK(se.n_modes == 6);
  CX_CHECK(!se.m_reduced.empty());

  // Solve the small reduced generalized eigenproblem K_r x = λ M_r x (dense) and compare
  // the lowest frequencies to the full model. Build LinearSystems from the reduced dense
  // matrices and reuse the eigensolution engine.
  const std::size_t d = se.dim();
  const fem::LinearSystem Kr = system_from_dense(se.k_reduced, d);
  const fem::LinearSystem Mr = system_from_dense(se.m_reduced, d);
  const numerics::EigenBasis red = numerics::extract_modes(Kr, Mr, 6);

  // The lowest 3 reduced frequencies must approximate the full model's within a few
  // percent (Craig-Bampton is exact in the limit of all interior modes; with 6 modes the
  // low ones are captured tightly). Compare relative error.
  for (std::size_t i = 0; i < 3; ++i) {
    const Real ff = full.modes[i].frequency;
    const Real fr = red.modes[i].frequency;
    CX_CHECK(ff > 0.0);
    const Real rel = std::fabs(fr - ff) / ff;
    CX_CHECK(rel < 0.02);  // < 2% — reduced low modes track the full model
  }
}

// Parse the stock-CalculiX substructure.inp and compare the exported reduced stiffness
// to its .dat.ref *MATRIX TYPE=STIFFNESS block (Guyan static reduction, 60 retained
// DOFs). Skipped when the CalculiX test tree is absent.
void test_reference_substructure_deck() {
  const std::string base = "/home/leonardo/work/CalculiX/test/substructure";
  std::ifstream inp(base + ".inp");
  std::ifstream ref(base + ".dat.ref");
  if (!inp.good() || !ref.good()) {
    std::printf("[test_substructure] reference deck absent — skipping\n");
    return;
  }
  const Model m = cxpp::io::parse_inp_file(base + ".inp");
  CX_CHECK(m.procedure == Procedure::Substructure);
  CX_CHECK(m.retained_dofs.size() == 60);  // N1 = 20 nodes × 3 DOFs
  const numerics::Superelement se = numerics::generate_substructure(m);
  CX_CHECK(se.dim() == 60);
  CX_CHECK(se.m_reduced.empty());  // STIFFNESS-only output -> no mass

  // Read the reference lower-triangular matrix values from the .dat.ref MATRIX block.
  std::string line;
  bool in_matrix = false;
  std::vector<Real> ref_vals;
  while (std::getline(ref, line)) {
    if (line.find("*MATRIX TYPE=STIFFNESS") != std::string::npos) {
      in_matrix = true;
      continue;
    }
    if (!in_matrix) continue;
    std::istringstream ls(line);
    Real v;
    while (ls >> v) ref_vals.push_back(v);
  }
  // 60×61/2 = 1830 lower-triangular entries.
  CX_CHECK(ref_vals.size() == 60 * 61 / 2);

  // Compare our lower triangle to the reference (relative tolerance on each entry; the
  // matrix spans ~1e2..1e6, so use a scaled tolerance against the reference magnitude
  // with an absolute floor for the near-zero entries).
  std::size_t idx = 0;
  Real max_rel = 0.0;
  for (std::size_t i = 0; i < 60; ++i)
    for (std::size_t j = 0; j <= i; ++j) {
      const Real ours = se.k_reduced[i * 60 + j];
      const Real theirs = ref_vals[idx++];
      const Real denom = std::max(std::fabs(theirs), 1.0);
      max_rel = std::max(max_rel, std::fabs(ours - theirs) / denom);
    }
  std::printf("[test_substructure] reference max relative entry error = %.3e\n",
              max_rel);
  CX_CHECK(max_rel < 1e-6);
}

}  // namespace

int main() {
  test_guyan_3dof_analytical();
  test_no_interior_is_identity();
  test_empty_retained_throws();
  test_export_format_lower_triangle();
  test_craig_bampton_approximates_full();
  test_reference_substructure_deck();
  CX_MAIN_RETURN();
}
