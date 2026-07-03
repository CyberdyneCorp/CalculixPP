#pragma once
#include <array>
#include <vector>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/types.hpp"

namespace cxpp::fem {

// Local element node indices (0-based) of a tetrahedral face (1..4), CalculiX
// ifacet ordering: 3 corners for C3D4, 3 corners + 3 midsides for C3D10.
std::vector<int> face_nodes(ElementType type, int face);

// Surface integrals over one element face, the shared kernel for face-based
// thermal boundary conditions (*FILM, *RADIATE). For a face with `n` nodes:
//   `gnode[i]`    global (mesh) node index of face node i
//   `load[i]`     = ∫ N_i dA               (consistent nodal "1" over the face)
//   `mass[i*n+j]` = ∫ N_i N_j dA           (face mass matrix, row-major)
// Film adds h*mass to the conduction matrix and h*T_sink*load to the rhs; radiation
// linearizes its q(T) into the same load/mass shapes. (spec: heat-transfer.)
struct FaceSurface {
  std::vector<Index> gnode;  // global node indices, size n
  std::vector<Real> load;    // ∫N_i dA, size n
  std::vector<Real> mass;    // ∫N_i N_j dA, size n*n row-major
};
FaceSurface face_surface_integrals(const Mesh& mesh, Index elem_index, int face);

// One face quadrature point in PHYSICAL space, for geometry that needs the surface
// sampled (cavity-radiation view factors): `x` is the point, `n` the UNIT outward
// normal (element -> outside, same sense as the *DLOAD pressure/*DFLUX face normal),
// and `w` the area weight so Σ w = face area and Σ w*x / area = the area centroid.
struct FacePoint {
  Vec3 x{0, 0, 0};
  Vec3 n{0, 0, 0};
  Real w{0.0};
};

// The face quadrature points of one element face in physical space (position, unit
// outward normal, area weight). Reuses the same face topology/shape/Gauss machinery
// as the surface-integral loads. `sub` (>=1) tiles the face into sub x sub subcells and
// applies the base rule in each, refining the sampling (Σ w stays the face area); the
// cavity-radiation view-factor double-area quadrature uses sub>1 for accuracy on far/
// opposed patches, where a single 2x2 rule underestimates the 1/r^2 kernel.
std::vector<FacePoint> face_gauss_points(const Mesh& mesh, Index elem_index, int face,
                                         int sub = 1);

// One element face sampled at its natural (parametric) coordinates: the physical
// position and the two covariant surface tangents t_xi = dx/dxi, t_eta = dx/deta at
// that parametric point, using the SAME face topology/shape machinery as the surface
// loads. Triangular faces (tet faces, wedge end faces) use barycentric (xi,eta) on the
// unit reference triangle (l1 = 1-xi-eta); quad faces (hex faces, wedge side faces) use
// (xi,eta) on [-1,1]^2. The outward normal is t_xi x t_eta (same sense as the *DLOAD
// pressure / *DFLUX face normal). This is the reusable geometric primitive the contact
// node-to-surface projection iterates on: it drives the residual dx·t_a of the
// closest-point problem and, at the solution, yields the projection point and the
// local normal/tangent frame. (spec: contact-search — node-to-surface projection.)
struct FaceFrame {
  Vec3 x{0, 0, 0};     // physical position of the parametric point
  Vec3 t_xi{0, 0, 0};  // dx/dxi
  Vec3 t_eta{0, 0, 0}; // dx/deta
};
FaceFrame face_frame_at(const Mesh& mesh, Index elem_index, int face, Real xi,
                        Real eta);

// One element face evaluated at parametric (xi,eta) in the DEFORMED configuration
// (reference node coords + the nodal displacement `u`, size mesh.num_nodes()). Returns
// the face node global indices, the shape values N_i at (xi,eta), and the deformed
// frame (position + covariant tangents). This is the primitive the node-to-surface
// contact operator uses each Newton iteration: it projects the deformed slave node onto
// the deformed master face and distributes the contact reaction to the master face nodes
// via N_i. With u all-zero it reduces to face_frame_at on the reference mesh. (spec:
// contact-search — node-to-surface projection in the current configuration.)
struct FaceEval {
  std::vector<Index> gnode;  // face node global (mesh) indices, size nf
  std::vector<Real> N;       // shape values N_i at (xi,eta), size nf
  Vec3 x{0, 0, 0};           // deformed physical position
  Vec3 t_xi{0, 0, 0};        // dx/dxi (deformed)
  Vec3 t_eta{0, 0, 0};       // dx/deta (deformed)
};
FaceEval face_eval_deformed(const Mesh& mesh, Index elem_index, int face, Real xi,
                            Real eta, const std::vector<Vec3>& u);

// Global thermal load vector (size num_nodes) at step fraction `lambda` in [0,1]:
// concentrated nodal heat flux (*CFLUX) plus consistent nodal fluxes from
// distributed surface flux (*DFLUX S<face>), each scaled by amplitude_factor.
// (spec: heat-transfer-analysis — *CFLUX / *DFLUX into the thermal rhs.)
std::vector<Real> thermal_load_vector(const Model& model, Real lambda);
std::vector<Real> thermal_load_vector(const Model& model);

// Global thermal-strain equivalent load vector (size 3*num_nodes) for a mechanical
// step with thermal expansion (spec: heat-transfer — thermal expansion coupling).
// For each element with *EXPANSION(alpha,Tref) and the model's applied_temperature
// field, adds f_th = ∫ Bᵀ D (alpha (T - Tref)) dV, so heating a constrained body
// induces thermal stress and a free body expands stress-free. Elements without
// expansion or without an applied temperature contribute nothing; the returned
// vector is all-zero when the model has no thermal strain (mechanical path
// unchanged). `lambda` scales the temperature change by the step fraction so an
// unamplified thermal load ramps 0->1 like the mechanical loads.
std::vector<Real> thermal_strain_load_vector(const Model& model, Real lambda);
std::vector<Real> thermal_strain_load_vector(const Model& model);

// Global external load vector (size 3*num_nodes) at full magnitude: concentrated
// loads (*CLOAD), consistent nodal forces from pressure (*DLOAD P<face> / *DSLOAD),
// and body loads (*DLOAD GRAV / CENTRIF). Amplitude references are evaluated at the
// end of the step (factor 1), so this is the fully-applied load.
std::vector<Real> external_load_vector(const Model& model);

// Same, but scaled per load by the step fraction `lambda` in [0,1]: each load is
// multiplied by amplitude_factor(load.amplitude, lambda * total), which reduces to
// `lambda` for loads without an amplitude. Used by the nonlinear driver so
// time-varying (*AMPLITUDE) loads are applied at the correct step fraction, while
// unamplified loads keep the linear 0->1 ramp. (spec: nonlinear-solution-control —
// Amplitude-driven time stepping.)
std::vector<Real> external_load_vector(const Model& model, Real lambda);

}  // namespace cxpp::fem
