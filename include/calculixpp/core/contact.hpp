#pragma once
#include <string>

#include "calculixpp/core/types.hpp"

// User-facing contact modeling cards (spec: contact — *CONTACT PAIR, *SURFACE
// INTERACTION, *SURFACE BEHAVIOR). These are the parsed, stored representation; the
// geometric search + penalty operator that consumes them lives in fem/contact.{hpp,cpp}
// and is assembled into the Newton tangent/residual by the nonlinear driver. A deck with
// no *CONTACT PAIR carries none of these and is byte-for-byte the pre-contact path.
namespace cxpp {

// Normal pressure-overclosure law from *SURFACE BEHAVIOR. The contact pressure p is a
// function of the signed normal gap g (g < 0 = penetration, g > 0 = clearance):
//   Hard:        p = -k * g for g < 0, else 0. k is a large penalty stiffness (~c*E),
//                so overclosure is (near) zero. `k` and (optional) `c0`/`p0` unused.
//   Linear:      p = -k * g for g < 0, else 0. Same form as hard but with a USER slope
//                k (the *SURFACE BEHAVIOR, PRESSURE-OVERCLOSURE=LINEAR slope), a softer
//                spring than the hard default. (k = slope.)
//   Exponential: p = p0 * (exp(c * (-g/c0)) - 1) / (exp(c) - 1) style softening. We use
//                the CalculiX two-parameter form p(g) = p0/(exp(1)-1) *
//                (exp(-(c0? )...)) — see the law in fem/contact.cpp. Parameters are
//                (p0 = pressure at zero clearance, c0 = clearance at which p ~ 0).
// `k` defaults to 0 meaning "auto" for the hard law (the operator picks c*E from the
// adjacent element stiffness); linear/exponential always carry explicit parameters.
struct SurfaceBehavior {
  enum class Law { Hard, Linear, Exponential };
  Law law{Law::Hard};
  Real k{0.0};   // Hard/Linear penalty slope (0 = auto for hard)
  Real p0{0.0};  // Exponential: pressure at zero clearance (c0 in the CalculiX card)
  Real c0{0.0};  // Exponential: clearance at which the pressure has ~vanished
};

// Coulomb friction from *FRICTION on a *SURFACE INTERACTION (spec: contact — tangential
// contact with stick/slip). The tangential traction sticks (elastic stick spring of
// stiffness `stick`) until it reaches the Coulomb limit mu * p (p the normal contact
// pressure force), then slips (return-mapping onto the friction cone: the tangential
// traction is capped at mu * p and points along the tangential slip direction). `stick`
// is the *FRICTION stick (lambda) stiffness (force per unit tangential slip); 0 -> the
// operator auto-sizes it from the normal penalty. `has` distinguishes an interaction with
// no *FRICTION (frictionless normal contact) from mu == 0 (given but zero).
struct Friction {
  Real mu{0.0};     // Coulomb friction coefficient
  Real stick{0.0};  // tangential stick stiffness (0 = auto from the normal penalty)
  bool has{false};  // true once a *FRICTION card set it
};

// Thermal contact conductance from *GAP CONDUCTANCE on a *SURFACE INTERACTION (spec:
// contact — thermal contact conductance and gap heat generation). Heat crosses a closed
// contact interface by conductance: the interface flux per unit area is q = h_c (T_slave
// - T_master), so heat flows from the hotter to the colder surface. `h` is the gap
// conductance coefficient (heat per area per unit temperature drop). `has` distinguishes
// an interaction with no *GAP CONDUCTANCE (thermally open — no cross-interface heat) from
// h == 0. Pressure/clearance dependence of h_c is a later refinement; this slice uses the
// constant coefficient evaluated wherever the surfaces are in contact (gap <= adjust).
struct GapConductance {
  Real h{0.0};      // gap conductance coefficient (per unit area per unit temperature)
  bool has{false};  // true once a *GAP CONDUCTANCE card set it
};

// Heat generated in the contact gap from *GAP HEAT GENERATION on a *SURFACE INTERACTION
// (spec: contact — gap heat generation). A heat flux per unit contacting area is generated
// at the interface (e.g. frictional/electrical) and split between the two surfaces. `q` is
// the generated flux per unit area (total; half deposited on each surface). `has`
// distinguishes an interaction with no card from q == 0.
struct GapHeatGeneration {
  Real q{0.0};      // generated heat flux per unit contacting area (total)
  bool has{false};  // true once a *GAP HEAT GENERATION card set it
};

// A *SURFACE INTERACTION: a named container the *CONTACT PAIR references, holding the
// normal surface behavior and (optional) Coulomb friction. `has_behavior` distinguishes
// an interaction that never got a *SURFACE BEHAVIOR (defaults to hard). `conductance` /
// `heat_generation` are the optional thermal-contact cards (inert on a purely mechanical
// or a thermally-open interface).
struct SurfaceInteraction {
  std::string name;
  SurfaceBehavior behavior{};
  bool has_behavior{false};
  Friction friction{};
  GapConductance conductance{};
  GapHeatGeneration heat_generation{};
};

// A *CONTACT PAIR (spec: contact — contact pairs). Node-to-surface is the implemented
// formulation: the slave *SURFACE (TYPE=NODE, or the nodes of an element surface) is the
// set of slave nodes projected onto the master element-face *SURFACE. `interaction` names
// the *SURFACE INTERACTION whose behavior supplies the pressure-overclosure law.
struct ContactPair {
  enum class Formulation { NodeToSurface, SurfaceToSurface };
  Formulation formulation{Formulation::NodeToSurface};
  std::string slave_surface;
  std::string master_surface;
  std::string interaction;
  // Search/adjust distance for the proximity grid (ADJUST=/small-sliding search). 0 ->
  // the operator picks a default from the master-face size.
  Real search{0.0};
  // Initial clearance from *CLEARANCE (spec: contact — contact modifiers). A uniform
  // initial gap ADDED to the geometric signed gap, so the contact operator sees
  // g_eff = g_geometric + clearance regardless of the as-meshed geometry: a positive
  // value opens the interface (needs that much closure before contact), a negative one
  // pre-loads it. 0 (no *CLEARANCE) leaves the geometric gap unchanged. `has_clearance`
  // distinguishes an explicit clearance of 0 from none (both behave identically here).
  Real clearance{0.0};
  bool has_clearance{false};
};

}  // namespace cxpp
