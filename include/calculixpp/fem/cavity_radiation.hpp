#pragma once
#include <vector>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/types.hpp"
#include "calculixpp/fem/loads.hpp"

// Gray-body cavity (surface-to-surface) radiation (spec: heat-transfer-analysis —
// radiation with view factors). Every *RADIATE ...,CR face is a PATCH of a diffuse-
// gray enclosure. Unlike surface-to-ambient radiation (which sees a fixed sink), a
// patch exchanges heat with every OTHER patch through the geometric view factors
// F_ij, weighted by the emissivities. This module computes the view factors and the
// per-patch net radiative heat flux (plus its temperature tangent) so the thermal
// Newton driver can add the coupling to the system, exactly parallel to the existing
// surface-to-ambient add_radiation.
namespace cxpp::fem {

// One radiating cavity patch: an element face reduced to the geometric quantities the
// view-factor kernel and the gray-body assembly need. `centroid`, `area` and the unit
// outward `normal` drive the double-area view-factor quadrature; `surf` (the shared
// face load/mass integrals) scatters the patch's net flux back onto its nodes; `emis`
// is the patch emissivity. `gp_x`/`gp_n`/`gp_w` are the face Gauss points in physical
// space (position, outward normal, area weight) used for the F_ij quadrature.
struct CavityPatch {
  Vec3 centroid{0, 0, 0};
  Vec3 normal{0, 0, 0};   // unit outward normal (element -> cavity)
  Real area{0.0};
  Real emis{0.0};
  FaceSurface surf;               // gnode / load (∫N) / mass (∫N N) over the face
  std::vector<Vec3> gp_x;         // Gauss-point positions
  std::vector<Vec3> gp_n;         // Gauss-point unit outward normals
  std::vector<Real> gp_w;         // Gauss-point area weights (Σ = patch area)
};

// The assembled cavity: the patches and the row-normalized view-factor matrix F
// (row-major, size n*n, F[i*n+j] = F_ij, the fraction of energy leaving patch i that
// reaches patch j). Reciprocity A_i F_ij = A_j F_ji holds to quadrature accuracy; each
// row is rescaled so Σ_j F_ij = 1 (a closed enclosure — the summation rule).
struct Cavity {
  std::vector<CavityPatch> patches;
  std::vector<Real> F;  // n*n row-major view-factor matrix, rows sum to 1
  int n{0};
};

// Build the cavity from the model's *RADIATE ...,CR faces: gather each patch's
// geometry + shared face integrals and compute the view-factor matrix by direct
// double-area quadrature F_ij = (1/A_i) ∫_Ai ∫_Aj cos_i cos_j / (pi r^2) dA_j dA_i,
// then row-normalize to Σ_j F_ij = 1 (closed enclosure). Self-view (i==j) is skipped
// (planar/convex patches see no part of themselves) and folded into the normalization.
// Returns an empty cavity (n == 0) when the model has no cavity radiation.
Cavity build_cavity(const Model& model);

// The per-patch net radiative heat flow Q_i (energy per unit time LEAVING patch i into
// the enclosure, a heat LOSS) and its tangent dQ_i/dT_j at the current patch mean
// absolute temperatures `tabs` (size n). Solves the gray-body radiosity system
//   (I - (1-eps) F) J = eps sigma T^4,   Q_i = A_i eps_i/(1-eps_i) (sigma T_i^4 - J_i),
// (with the black-body limit eps -> 1 handled directly, Q_i = A_i (sigma T_i^4 -
// Σ_j F_ij sigma T_j^4)). `sigma` is Stefan-Boltzmann. Fills `Q` (size n) and `dQdT`
// (row-major n*n, dQ_i/dT_j). Used by the thermal driver to build the Newton update.
void cavity_heat_flow(const Cavity& cav, const std::vector<Real>& tabs, Real sigma,
                      std::vector<Real>& Q, std::vector<Real>& dQdT);

}  // namespace cxpp::fem
