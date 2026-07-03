#pragma once
#include <vector>

#include "calculixpp/core/contact.hpp"
#include "calculixpp/core/model.hpp"
#include "calculixpp/core/results.hpp"
#include "calculixpp/core/types.hpp"
#include "calculixpp/fem/assembly.hpp"
#include "calculixpp/fem/contact_search.hpp"

// Penalty node-to-surface contact operator (spec: contact — contact pairs / surface
// behavior; contact-search tasks 1.4/1.5). This is where the geometric search
// (contact_search.hpp) meets the pressure-overclosure law (SurfaceBehavior) and the
// Newton spine: for each active slave-node/master-face pair with penetration g < 0 it
// produces a normal contact force (into the residual) and its consistent penalty tangent
// (into the reduced free/free system), so contact is just another contribution to the
// Newton tangent + residual (spec: contact-search — contact operator generation through
// the ComputeBackend). A deck with no *CONTACT PAIR never builds any of this.
//
// Math reference (reimplemented, not copied): CalculiX gencontelem_n2f.f (n2f penalty
// spring) and surfacebehaviors.f (the overclosure laws) — algorithm reference only.
namespace cxpp::fem {

// One resolved contact pair ready to evaluate each Newton iteration: the slave node
// indices (dense mesh node indices), the master faces, a spatial grid over them, and the
// normal behavior law. Built once per solve by build_contact_pairs (the geometry the
// grid indexes is the REFERENCE mesh; small-sliding node-to-surface re-projects the
// deformed slave node against the reference master faces plus the current master motion).
struct ResolvedContactPair {
  std::vector<Index> slave_nodes;      // dense mesh node indices of the slave surface
  std::vector<MasterFace> master_faces;
  MasterGrid grid;
  SurfaceBehavior behavior;
  Friction friction;  // Coulomb friction (frictionless when friction.has == false)
  Real stick{0.0};    // effective tangential stick stiffness (auto or from *FRICTION)
  Real penalty{0.0};  // effective nodal penalty stiffness kappa (auto or from behavior)
  // Initial clearance from *CLEARANCE, added to the geometric signed gap the operator
  // sees (g_eff = g_geometric + clearance). 0 when the pair has no *CLEARANCE. (spec:
  // contact — contact modifiers / *CLEARANCE.)
  Real clearance{0.0};
  // Adjust band: a small positive normal-gap tolerance. A slave node within [0, adjust)
  // of the master face (near-contact) contributes its penalty STIFFNESS to the tangent so
  // the interface has no spurious rigid-body mode at first touch, while the contact FORCE
  // is still zero until the node actually penetrates (gap < 0). Auto-sized from the master
  // element size. (spec: contact — hard contact with near-zero allowed overclosure.)
  Real adjust{0.0};
};

// Resolve every *CONTACT PAIR of the model into a ResolvedContactPair (slave nodes +
// master faces + grid + behavior + penalty stiffness). Node-to-surface only; a
// surface-to-surface pair throws (mortar is a later workstream). The penalty stiffness is
// taken from the behavior (Linear/Hard k, or an auto estimate ~ c * E * h from the
// adjacent element for a hard law with k == 0). (spec: contact — pairs consume named
// surfaces + a surface interaction.)
std::vector<ResolvedContactPair> build_contact_pairs(const Model& model);

// Evaluate the contact operator at the current full nodal displacement `u` (size
// mesh.num_nodes()) and add its contributions to the ALREADY-REDUCED Newton system:
//   - the normal contact force is accumulated into `f_int` (full nodal DOF length,
//     size 3*num_nodes) as ∂Π/∂u so the driver's residual r = f_ext - f_int drives the
//     slave out of penetration;
//   - the penalty tangent k (n⊗n across the slave↔master DOFs) is scattered into the
//     reduced COO system `sys` (rows/cols/vals) through its constraint transform, with
//     any prescribed-DOF columns moved to sys.rhs — exactly the congruence path the
//     element tangents use.
// When the pair carries Coulomb *FRICTION (friction.has, mu > 0), a penetrating slave node
// also gets a TANGENTIAL traction: the elastic stick predictor -stick*g_t (g_t the slave's
// tangential slip vs its master footprint) return-mapped onto the friction cone |t| <= mu*p
// (stick inside, slip capped at mu*p outside), added to f_int with its consistent regime
// tangent — the stick/slip state is re-decided each call (active-set within the Newton loop).
// Only penetrating pairs (g < 0) contribute (the active set); a released (tensile) pair
// contributes nothing. Returns the number of active slave nodes (diagnostics). The
// reference geometry of the master faces plus the current displacement `u` give the
// deformed positions used for projection. (spec: contact-search — contact contributes to
// tangent and residual; contact — normal + tangential contact.)
int add_contact(const Model& model, const std::vector<ResolvedContactPair>& pairs,
                const std::vector<Vec3>& u, LinearSystem& sys, std::vector<Real>& f_int);

// One resolved thermal-contact pair (spec: contact — thermal contact conductance / gap
// heat generation). Shares the SAME node-to-surface search as mechanical contact (task
// 1.6 — one interface geometry): the slave nodes, the master faces + grid, and the
// conductance/heat-generation coefficients. `slave_area[i]` is the tributary contact area
// of slave_nodes[i] (∫N dA over the slave surface, or master-area/count when the slave is
// a node surface) so the interface conductance/heat is area-consistent.
struct ResolvedThermalContactPair {
  std::vector<Index> slave_nodes;
  std::vector<Real> slave_area;  // tributary area per slave node (aligned with slave_nodes)
  std::vector<MasterFace> master_faces;
  MasterGrid grid;
  Real h{0.0};        // gap conductance coefficient (per area per temperature)
  Real q_gap{0.0};    // generated gap heat flux per unit area (total)
  Real adjust{0.0};   // in-contact gap band (same near-contact tolerance as mechanical)
};

// Resolve every *CONTACT PAIR whose interaction carries *GAP CONDUCTANCE or *GAP HEAT
// GENERATION into a ResolvedThermalContactPair. Empty when the model has no thermal
// contact, so a plain thermal deck is unchanged. (spec: contact — thermal contact.)
std::vector<ResolvedThermalContactPair> build_thermal_contact_pairs(const Model& model);

// Add the thermal-contact conductance + gap-heat-generation contributions to a scalar
// thermal system (spec: contact — thermal contact conductance and gap heat generation).
// For each slave node in contact (projected gap <= adjust, in the geometry given by the
// nodal displacement `u` — pass an all-zero field for a pure thermal deck) it adds the
// conductance coupling  Q = h A_s (T_slave - Σ N_i T_master_i)  to the operator: a
// symmetric conductance block over the slave node and the master face nodes (added to
// `k_rows/k_cols/k_vals` COO over node indices), and, for gap heat generation, the source
// q_gap A_s / 2 deposited on the slave node and split over the master face nodes (added to
// `source`, size num_nodes). The temperatures cancel in the conductance block so it is
// heat that flows h A (T_s - T_m); no temperature field is needed here (the matrix carries
// the coupling). Returns the number of slave nodes found in contact. `u` is size
// mesh.num_nodes(); `source` is size mesh.num_nodes(). (spec: contact — gap state drives
// thermal conductance; shared search for thermal contact.)
int add_thermal_contact(const Model& model,
                        const std::vector<ResolvedThermalContactPair>& pairs,
                        const std::vector<Vec3>& u, std::vector<Index>& k_rows,
                        std::vector<Index>& k_cols, std::vector<Real>& k_vals,
                        std::vector<Real>& source);

// Recover the per-slave-node contact result (CSTR) at the converged displacement `u`
// (spec: contact — contact output). For each slave node it re-runs the node-to-surface
// projection in the deformed config, and reports the status (CLOSED when penetrating),
// the normal contact PRESSURE p = F_normal / tributary_area (>= 0), the signed gap, and
// the tangential friction traction magnitude. Uses the same penalty/friction laws as the
// operator (add_contact), so the reported pressure is the assembled contact pressure. The
// tributary area is the master-face-based slave area (shared with thermal contact). Empty
// when the model has no contact pair. (spec: contact — CSTR to .dat/.frd.)
std::vector<ContactPoint> recover_contact(const Model& model,
                                          const std::vector<ResolvedContactPair>& pairs,
                                          const std::vector<Vec3>& u);

// The normal contact pressure-force weight for a signed gap `g` under `behavior` with
// effective penalty stiffness `kappa`: the magnitude of the restoring force per active
// slave node (>= 0), zero for a positive clearance (g >= 0). Exposed for unit testing of
// the overclosure laws. (spec: contact-search — normal pressure from the overclosure
// law.)
Real contact_force_weight(const SurfaceBehavior& behavior, Real kappa, Real g);

// Coulomb friction return-mapping of a tangential-traction TRIAL vector (spec: contact —
// stick/slip). Given the elastic-predictor tangential traction `t_trial` (= -stick * slip,
// already projected into the tangent plane), the friction coefficient `mu` and the normal
// contact pressure force `p_normal` (>= 0), returns the admissible tangential traction:
//   - |t_trial| <= mu*p_normal  -> STICK: return t_trial unchanged; `slipping` = false.
//   - |t_trial| >  mu*p_normal  -> SLIP: return-map onto the friction cone, scaling
//     t_trial down to magnitude mu*p_normal (same direction); `slipping` = true.
// mu == 0 or p_normal == 0 caps the traction at 0 (pure slip / no normal pressure).
// Exposed for unit-testing the stick<->slip transition. (spec: contact — Coulomb cone.)
Vec3 friction_return_map(const Vec3& t_trial, Real mu, Real p_normal, bool& slipping);

}  // namespace cxpp::fem
