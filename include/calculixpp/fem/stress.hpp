#pragma once
#include <vector>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/results.hpp"
#include "calculixpp/core/types.hpp"
#include "calculixpp/fem/assembly.hpp"

namespace cxpp::fem {

// Global internal force vector f_int = sum_e Ke_e u_e (size 3*num_nodes, DOF order
// [u,v,w] per node). For linear elasticity Ke_e is the element stiffness, so this is
// the exact internal force that the Newton residual r = lambda*f_ext - f_int uses;
// it is structured to accept a nonlinear tangent/stress later. `u` is the current
// nodal displacement (all nodes). (spec: nonlinear-solution-control — residual.)
std::vector<Real> internal_force(const Model& model, const std::vector<Vec3>& u);

// Given a solved displacement field (fields.displacement, all nodes), recover the
// averaged nodal stresses and the reaction forces. Integration-point stresses are
// extrapolated to element nodes (C3D4: constant; C3D10: a4 corner extrapolation +
// midside averaging) and averaged across elements. RF = f_int - f_ext.
// (spec: static-analysis / results-output — reimplemented from CalculiX
// resultsmech.f / extrapolate.f, not copied.)
void recover_fields(const Model& model, StaticFields& fields);

// Per-element per-integration-point reference stress (Voigt6 {xx,yy,zz,xy,xz,yz})
// for a supplied full nodal displacement `u` (spec: geometric-stiffness — two-step
// prestress driver). `out[e][q]` is the Gauss-point Cauchy stress of element e used
// to feed assemble_geometric_stiffness. It shares the exact integration-point math as
// recover_fields (same Gauss rule, same thermal-corrected σ = D (eps_mech - eps_th)),
// factored out so the buckling prestress recovery reuses the shipped path; recover_fields
// itself is left byte-identical. Deactivated elements (*MODEL CHANGE, REMOVE) get an
// empty per-element vector. (ref: resultsini.f / resultsmech.f.)
std::vector<std::vector<Voigt6>> recover_gauss_stress(const Model& model,
                                                      const std::vector<Vec3>& u);

// Plasticity-aware field recovery. Identical to recover_fields for elastic elements,
// but the integration-point stress/strain and the internal force come from the
// committed material-point state (radial-return stress, not linear D*strain), so a
// plastic solve reports the true stress on the hardening curve. `mp` must be the
// driver's converged MaterialPoints (committed history), aligned with mesh.elements().
// (spec: material-models / results-output — stress from the committed plastic state.)
void recover_fields(const Model& model, StaticFields& fields, MaterialPoints& mp);

}  // namespace cxpp::fem
