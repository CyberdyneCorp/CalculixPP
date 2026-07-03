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
