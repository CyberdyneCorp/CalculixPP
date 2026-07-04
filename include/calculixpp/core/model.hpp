#pragma once
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "calculixpp/core/amplitude.hpp"
#include "calculixpp/core/connector.hpp"
#include "calculixpp/core/constraint.hpp"
#include "calculixpp/core/constraint_cards.hpp"
#include "calculixpp/core/contact.hpp"
#include "calculixpp/core/material.hpp"
#include "calculixpp/core/mesh.hpp"
#include "calculixpp/core/section.hpp"
#include "calculixpp/core/solution_control.hpp"
#include "calculixpp/core/submodel.hpp"

namespace cxpp {

// Single-point constraint: prescribed nodal DOF (spec: loads-and-boundary-conditions).
// `amplitude` (optional *AMPLITUDE name) scales the prescribed value over step time.
struct Spc {
  Index node_id{};
  int comp{};      // 1..3 (x,y,z)
  Real value{0.0};
  std::string amplitude;  // empty -> default linear ramp
  // *BOUNDARY, FIXED: hold the DOF at its currently-attained (deformed) value going
  // into the step, rather than at `value`. Only meaningful in multi-step analysis —
  // the step-loop driver (numerics::solve_multistep_static) substitutes the carried
  // displacement of the DOF for `value` at solve time. In a single-step solve there is
  // no prior deformation so FIXED reduces to value 0 (the DOF's start state), i.e. the
  // single-step path is unchanged. (spec: multi-step analysis — FIXED boundary.)
  bool fixed{false};
  // *BOUNDARY, SUBMODEL: this DOF is DRIVEN from a global solution (spec: submodeling).
  // `value` is a placeholder — the submodel driver overwrites it with the global
  // displacement interpolated at the node's location before the local solve. A non-driven
  // SPC (the default) is a literal prescribed value, so the ordinary path is unchanged.
  // Declared AFTER `fixed` so the existing `Spc{...,fixed}` aggregate inits stay valid.
  bool driven{false};
};

// Concentrated load on a nodal DOF (*CLOAD). `amplitude` scales the value over
// step time (empty -> default linear ramp).
struct Cload {
  Index node_id{};
  int comp{};      // 1..3
  Real value{0.0};
  std::string amplitude;
};

// Distributed pressure on an element face (*DLOAD P<face> / *DSLOAD). Positive
// pressure acts into the element (opposite the outward face normal). Face 1..4.
// `amplitude` scales the pressure over step time (empty -> default linear ramp).
struct Dload {
  Index elem_id{};
  int face{};      // 1..4
  Real pressure{0.0};
  std::string amplitude;
};

// Body load applied over element volume (*DLOAD label GRAV or CENTRIF), as
// consistent nodal forces integrated from rho * N over the element (needs
// *DENSITY). GRAV: `magnitude` is g, `a`..`c` the (normalized) direction.
// CENTRIF: `magnitude` is omega^2, `p` the axis point, `a`..`c` the axis
// direction. `amplitude` scales `magnitude` over step time.
struct BodyLoad {
  enum class Kind { Gravity, Centrifugal };
  Kind kind{Kind::Gravity};
  std::string elset;      // element set the load acts on (empty -> all elements)
  Real magnitude{0.0};    // GRAV: g;  CENTRIF: omega^2
  Vec3 dir{0, 0, 0};      // GRAV: gravity direction;  CENTRIF: axis direction
  Vec3 point{0, 0, 0};    // CENTRIF: a point on the rotation axis
  std::string amplitude;
};

// Prescribed nodal temperature for a heat-transfer step (spec:
// heat-transfer-analysis). In a heat step *BOUNDARY node,11[,11],value fixes the
// temperature DOF; *TEMPERATURE node,value does the same. `amplitude` scales the
// value over step time (empty -> default linear ramp), mirroring Spc.
struct TempBc {
  Index node_id{};
  Real value{0.0};
  std::string amplitude;
};

// Concentrated nodal heat flux (*CFLUX node,dof,value). The dof field (0 or 11 in
// CalculiX) selects the temperature DOF and is not stored; the flux adds directly
// to the thermal rhs. `amplitude` scales the value over step time.
struct Cflux {
  Index node_id{};
  Real value{0.0};
  std::string amplitude;
};

// Distributed surface heat flux on an element face (*DFLUX elem,S<face>,magnitude).
// Positive `flux` is heat flowing INTO the element across the face (integrated as
// ∫ N_a * flux dA over the face). Face 1..6 depending on element topology.
// `amplitude` scales the flux over step time. (spec: heat-transfer-analysis.)
struct Dflux {
  Index elem_id{};
  int face{};       // 1..6 (tet 1..4)
  Real flux{0.0};
  std::string amplitude;
};

// Convective film boundary condition on an element face (*FILM elem,F<face>,
// T_sink,h). The face loses heat q_film = h (T - T_sink) per unit area, which adds
// h ∫N_a N_b dA to the conduction matrix and h T_sink ∫N_a dA to the thermal rhs.
// `amplitude` scales the sink temperature over step time. Face 1..6 (tet 1..4).
// (spec: heat-transfer-analysis — convective film loads.)
struct Film {
  Index elem_id{};
  int face{};        // 1..6 (tet 1..4)
  Real sink_temp{0.0};
  Real h{0.0};       // film (convection) coefficient
  std::string amplitude;
};

// Radiation on an element face (*RADIATE). Two modes, distinguished by `cavity`:
//   surface-to-ambient (R<face>): the face radiates q = eps sigma (T^4 - T_amb^4)
//     per unit area against a fixed ambient/environment temperature;
//   cavity (R<face>CR): the face is a patch of a gray-body enclosure — the net
//     radiative flux depends on the temperatures AND emissivities of every other
//     cavity patch through the geometric view factors F_ij (see cavity_radiation).
// Both are in absolute temperature T_abs = T - absolute_zero, linearized within the
// thermal Newton iteration. `emissivity` is eps in [0,1]; the Stefan-Boltzmann
// constant and absolute-zero offset come from the model's PhysicalConstants.
// `amplitude` scales the ambient temperature over step time (unused for cavity).
// (spec: heat-transfer-analysis — radiation with view factors.)
struct Radiate {
  Index elem_id{};
  int face{};          // 1..6 (tet 1..4)
  Real ambient_temp{0.0};
  Real emissivity{0.0};
  bool cavity{false};  // true -> R<face>CR gray-body cavity patch (view factors)
  std::string amplitude;
};

// Which linear solver the step requested (via SOLVER= on *STATIC). Maps the
// CalculiX solver names onto the paths CalculiX++ implements. When SOLVER= is
// unspecified the request is Auto: the numerics layer picks sparse-direct for
// small/medium systems and IC0-preconditioned CG for large ones (memory-scalable
// on 3D/mobile). (spec: linear-algebra-and-solvers 9.2/9.3.)
enum class RequestedSolver {
  Auto,    // default; size-based direct-vs-CG choice (see resolve_solver_kind)
  Direct,  // SPOOLES/PARDISO/PASTIX/DIRECT map here (sparse spsolve)
  CG,      // ITERATIVE*/CG (IC0-preconditioned cg)
};

// Analysis procedure selected by the step's procedure card. Static is the default
// (mechanical *STATIC); HeatTransferSteady is a *HEAT TRANSFER, STEADY STATE step,
// HeatTransferTransient a transient *HEAT TRANSFER step (backward-Euler time
// integration of the capacitance term). Both are scalar temperature fields.
// (spec: heat-transfer-analysis.)
enum class Procedure {
  Static,
  HeatTransferSteady,
  HeatTransferTransient,
  // *COUPLED TEMPERATURE-DISPLACEMENT: temperature + displacement solved together.
  // The one-way (sequential) scheme solves the steady thermal field first, then
  // applies its thermal strain to a mechanical solve. (spec: heat-transfer — coupled.)
  Coupled,
  // *FREQUENCY: natural-frequency / mode-shape extraction — the generalized symmetric
  // eigenproblem K x = λ M x for the lowest N modes (spec: modal-and-buckling-analysis).
  Frequency,
  // *BUCKLE: linear-buckling load-factor extraction — the pencil (K + λ K_geo) φ = 0
  // where K_geo is the geometric (initial-stress) stiffness about a reference-load
  // prestress state. Returns the lowest positive load factors + buckling mode shapes
  // (spec: modal-and-buckling-analysis — *BUCKLE; geometric-stiffness).
  Buckling,
  // *MODAL DYNAMIC: transient response by modal superposition over a preceding
  // *FREQUENCY basis — each decoupled modal SDOF integrated, recombined to the physical
  // response (spec: dynamic-analysis — modal dynamic).
  ModalDynamic,
  // *STEADY STATE DYNAMICS: harmonic response amplitude/phase over a frequency sweep by
  // modal superposition with damping (spec: dynamic-analysis — steady-state dynamics).
  SteadyStateDynamics,
  // *COMPLEX FREQUENCY: damped complex modes by proportional (Rayleigh / modal) damping,
  // reduced onto the preceding *FREQUENCY basis — the reduced quadratic eigenproblem
  // (λ²M + λC + K)x = 0 projected onto Φ and solved as a small 2·nev companion pencil
  // (spec: modal-and-buckling-analysis — complex frequency). This is the option-(B)
  // proportional-damping path and is deliberately NOT the CalculiX CORIOLIS gyroscopic
  // problem; CORIOLIS/FLUTTER decks are parsed then rejected. (numerics/eigensolution.)
  ComplexFrequency,
  // *DYNAMIC: direct time integration of M a + C v + K u = f(t) by the implicit HHT-α
  // scheme in physical coordinates (spec: dynamic-analysis — direct-integration
  // dynamics). Carries inertia + nonlinear material/contact tangents together; the
  // undamped α=0 scheme conserves energy. (numerics/direct_dynamics.)
  Dynamic,
  // *SUBSTRUCTURE GENERATE: condense the model onto the retained (master) DOFs of
  // *RETAINED NODAL DOFS, producing a superelement. Static (Guyan) reduction of the
  // interior DOFs gives the reduced stiffness Schur complement; when a mass matrix is
  // present, Craig-Bampton adds fixed-interface normal modes and a reduced mass. The
  // reduced matrices are exported by *SUBSTRUCTURE MATRIX OUTPUT / *MATRIX ASSEMBLE.
  // (spec: substructure-generation.) (numerics/substructure.)
  Substructure,
};

// Solution scheme for a *COUPLED TEMPERATURE-DISPLACEMENT step (spec: heat-transfer
// — coupled). The two schemes solve the SAME coupled residual and agree to solver
// tolerance; they differ in how the thermal (T) and mechanical (u) blocks are
// combined.
//   - Monolithic: assemble and solve the 4-DOF/node (u,v,w,T) tangent in one Newton
//     system. The thermal->mechanical coupling enters through the thermal-strain
//     load; the mechanical->thermal block (K_Tu) is non-zero only when a genuine
//     two-way heat source exists (plastic-dissipation heating). With no two-way term
//     the system is block-triangular and the result equals the one-way solve exactly.
//   - Staggered: Gauss-Seidel — alternate a thermal solve and a mechanical solve with
//     an outer convergence check on the coupled fields. Converges in ONE pass when the
//     coupling is one-way (no mechanical->thermal feedback); iterates otherwise.
// The default is Staggered (its single-pass degeneracy is the historical one-way
// path). CalculiX's *COUPLED TEMPERATURE-DISPLACEMENT selects the monolithic tangent
// unless SOLUTIONS=STAGGERED; the parser maps that keyword.
enum class CoupledScheme { Staggered, Monolithic };

// A *MODEL CHANGE, TYPE=CONTACT PAIR record (spec: model-change — activate/deactivate
// a named contact pair between steps). Parsed and stored only in this thermal track;
// the contact assembly that consumes it lands in the contact workflow. `add` is true
// for ADD (activate), false for REMOVE (deactivate). The two surface names identify the
// pair (order as written on the *CONTACT PAIR card). (spec: model-change — contact pair.)
struct ContactPairChange {
  std::string surface_a;  // first *SURFACE name of the pair
  std::string surface_b;  // second *SURFACE name of the pair
  bool add{true};         // true = ADD (activate), false = REMOVE (deactivate)
};

// Physical constants for radiation and absolute-temperature conversion, from
// *PHYSICAL CONSTANTS (ABSOLUTE ZERO, STEFAN BOLTZMANN). Radiation works in
// absolute temperature T_abs = T - absolute_zero; sigma is the Stefan-Boltzmann
// constant. Zero sigma (no card) makes *RADIATE inert. (spec: heat-transfer.)
struct PhysicalConstants {
  Real absolute_zero{0.0};  // e.g. -273.15 for Celsius decks
  Real sigma{0.0};          // Stefan-Boltzmann constant
};

// Rayleigh (proportional) damping C = alpha*M + beta*K from *DAMPING, ALPHA=, BETA=
// (spec: dynamic-analysis — damping). Both zero (the default / no card) is undamped.
struct Rayleigh {
  Real alpha{0.0};  // mass-proportional coefficient
  Real beta{0.0};   // stiffness-proportional coefficient
  bool active{false};
};

// Base motion (support excitation) from *BASE MOTION (spec: dynamic-analysis — base
// motion). `dof` is the excitation translational direction (1=x,2=y,3=z); `magnitude`
// is the base acceleration amplitude applied along it. Inert unless active.
struct BaseMotion {
  int dof{0};             // excitation direction 1..3 (0 = none)
  Real magnitude{0.0};    // base acceleration amplitude
  bool active{false};
};

// The assembled analysis model for the linear-static slice.
class Model {
 public:
  Mesh mesh;
  std::unordered_map<std::string, Material> materials;
  std::vector<SolidSection> sections;
  std::vector<Spc> spcs;
  // *SUBMODEL cards (spec: submodeling). Non-empty routes a *STATIC deck through the
  // submodel driver, which fills the driven SPCs from a supplied global solution before
  // the local solve. Empty -> ordinary analysis (byte-for-byte the pre-submodel path).
  std::vector<SubmodelSpec> submodels;
  std::vector<Cload> cloads;
  std::vector<Dload> dloads;
  std::vector<BodyLoad> body_loads;
  std::unordered_map<std::string, Amplitude> amplitudes;

  // Thermal step data (spec: heat-transfer-analysis). Populated only for a
  // *HEAT TRANSFER step; inert on the mechanical path.
  Procedure procedure{Procedure::Static};
  std::vector<TempBc> temp_bcs;  // prescribed temperatures (*BOUNDARY dof 11 / *TEMPERATURE)
  std::vector<Cflux> cfluxes;    // concentrated nodal heat flux (*CFLUX)
  std::vector<Dflux> dfluxes;    // distributed surface heat flux (*DFLUX)
  std::vector<Film> films;       // convective film boundary (*FILM)
  std::vector<Radiate> radiates; // surface-to-ambient radiation (*RADIATE)
  PhysicalConstants physical;    // *PHYSICAL CONSTANTS (radiation sigma / abs zero)

  // Uniform initial temperature for a transient step (*INITIAL CONDITIONS,
  // TYPE=TEMPERATURE). Per-node overrides land in initial_temperature_by_node;
  // nodes absent from it fall back to initial_temperature. Steady state ignores it.
  Real initial_temperature{0.0};
  std::unordered_map<Index, Real> initial_temperature_by_node;

  // Applied nodal temperature field for a MECHANICAL step with thermal expansion
  // (spec: heat-transfer — coupled / thermal expansion). On a plain *STATIC deck the
  // parser fills this from *TEMPERATURE cards; a *COUPLED TEMPERATURE-DISPLACEMENT
  // step overwrites it with the solved thermal field before the mechanical solve.
  // Nodes absent from the map carry no temperature change (thermal strain zero). The
  // thermal strain of an element with *EXPANSION(alpha,Tref) is
  // eps_th = alpha (T - Tref) on the normal components. Empty -> pure mechanics.
  std::unordered_map<Index, Real> applied_temperature;

  // Solution scheme for a *COUPLED TEMPERATURE-DISPLACEMENT step (Monolithic vs.
  // Gauss-Seidel Staggered). Inert unless procedure == Coupled. (spec: heat-transfer
  // — coupled.)
  CoupledScheme coupled_scheme{CoupledScheme::Staggered};

  // Taylor-Quinney fraction beta in [0,1] for plastic-dissipation heating (the
  // mechanical->thermal two-way coupling term). A fraction beta of the incremental
  // plastic work density is deposited as a volumetric heat source in the coupled
  // thermal solve, so plastic straining warms the body, which changes the thermal
  // strain, which the coupled scheme iterates to convergence. Zero (the default)
  // makes the coupling one-way, so the coupled result equals the sequential solve
  // for pure thermal-stress problems. (spec: heat-transfer — coupled two-way.)
  Real taylor_quinney{0.0};

  // Per-element plastic-dissipation heat source Q_e (total watts deposited in element
  // e) for the coupled thermal solve, given the per-element accumulated equivalent
  // plastic strain `eqplastic_by_elem` (aligned with mesh.elements()) recovered from
  // a mechanical solve. Q_e = beta * sigma_y(ep_e) * ep_e * V_e, the Taylor-Quinney
  // fraction of the plastic work density integrated over the element volume. Returns
  // an all-zero vector when taylor_quinney == 0 or the model has no plasticity, so a
  // pure thermal-stress deck stays one-way. (spec: heat-transfer — coupled two-way.)
  std::vector<Real> plastic_dissipation_heat(
      const std::vector<Real>& eqplastic_by_elem) const;

  // Discrete/connector elements (spec: element-sections — *SPRING/*MASS/*DASHPOT).
  // Springs contribute to the static stiffness; masses and dashpots are stored for
  // later dynamics (Phase 4) and are inert in a static solve.
  std::vector<Spring> springs;
  std::vector<PointMass> point_masses;
  std::vector<Dashpot> dashpots;

  // Multi-point constraints (spec: constraints). `equations` are raw *EQUATION
  // linear relations (first term dependent); the higher-level cards below are
  // expanded into equations at assembly time by expand_constraints().
  std::vector<Equation> equations;
  std::vector<Mpc> mpcs;
  std::vector<RigidBody> rigid_bodies;
  std::vector<Coupling> couplings;
  std::vector<Tie> ties;

  // Element birth-death from *MODEL CHANGE, TYPE=ELEMENT (spec: model-change). Element
  // ids listed here are DEACTIVATED for the (single) step: they contribute no
  // stiffness, mass, internal force, body load, or thermal-strain load, and no stress.
  // A `TYPE=ELEMENT, ADD` for an id removes it from this set (reactivation). Because the
  // current model is single-step, reactivation is strain-free by construction — the
  // solve starts from the undeformed state, so a re-added element carries zero initial
  // strain. Cross-step birth-death (an element active in one step, removed in the next,
  // re-added strain-free relative to the deformed geometry) needs the multi-step engine
  // and is deferred (see openspec/BACKLOG.md / tasks.md 5.1). Empty -> every element
  // active (byte-for-byte the pre-model-change path).
  std::vector<Index> deactivated_elements;

  // *MODEL CHANGE, TYPE=CONTACT PAIR records, parsed and stored only here (consumed by
  // the contact workflow). (spec: model-change — contact pair activation.)
  std::vector<ContactPairChange> contact_pair_changes;

  // Contact modeling (spec: contact — *CONTACT PAIR / *SURFACE INTERACTION / *SURFACE
  // BEHAVIOR). `contact_pairs` are the parsed pairs; `surface_interactions` maps an
  // interaction name to its normal pressure-overclosure behavior. A deck with any
  // *CONTACT PAIR routes to the nonlinear (Newton) driver, where the penalty
  // node-to-surface operator (fem/contact.{hpp,cpp}) contributes to the tangent +
  // residual each iteration. Empty -> no contact (pre-contact path unchanged).
  std::vector<ContactPair> contact_pairs;
  std::unordered_map<std::string, SurfaceInteraction> surface_interactions;

  // True when the model carries at least one *CONTACT PAIR. Routes the analysis to
  // solve_nonlinear_static with the contact contribution; a contact-free deck keeps its
  // existing (linear or nonlinear-material) path unchanged. (spec: contact — a contact
  // deck is a nonlinear constraint problem.)
  bool has_contact() const { return !contact_pairs.empty(); }

  // Return a COPY of this model whose `contact_pairs` are filtered to the pairs currently
  // ACTIVE given the accumulated `*MODEL CHANGE, TYPE=CONTACT PAIR` records
  // (`contact_pair_changes`) — a pair is active unless a REMOVE (not later re-ADDed) names
  // its surface pair. A pair matches a change when its slave/master surfaces equal the
  // change's two surfaces in EITHER order (the *CONTACT PAIR data line and the *MODEL
  // CHANGE data line may list them in either order). With no contact-pair change every
  // pair stays active, so a deck without *MODEL CHANGE, TYPE=CONTACT PAIR is unchanged.
  // The multi-step engine calls this per step so an inactive pair adds no contact
  // search/force/constraint that step. (spec: model-change — activate/deactivate a
  // contact pair between steps.)
  Model with_active_contact_pairs() const;

  // The resolved pressure-overclosure behavior of a *CONTACT PAIR: the behavior of its
  // named *SURFACE INTERACTION, or a default hard behavior when the interaction carries
  // no *SURFACE BEHAVIOR (CalculiX defaults contact to hard). (spec: contact — surface
  // behavior normal.)
  SurfaceBehavior contact_behavior(const ContactPair& pair) const;

  // The resolved Coulomb friction of a *CONTACT PAIR: the friction of its named *SURFACE
  // INTERACTION, or a default frictionless Friction{} when the interaction carries no
  // *FRICTION card. (spec: contact — tangential contact / Coulomb friction.)
  Friction contact_friction(const ContactPair& pair) const;

  // The resolved thermal gap conductance / gap heat generation of a *CONTACT PAIR: the
  // corresponding card on its named *SURFACE INTERACTION, or a defaulted (has == false)
  // record when the interaction carries no such card. (spec: contact — thermal contact.)
  GapConductance contact_conductance(const ContactPair& pair) const;
  GapHeatGeneration contact_heat_generation(const ContactPair& pair) const;

  // True when the model has at least one *CONTACT PAIR whose interaction carries a
  // *GAP CONDUCTANCE or *GAP HEAT GENERATION card, so a thermal step must add the
  // cross-interface conductance/heat operator. A deck with no thermal-contact card keeps
  // its existing thermal path unchanged. (spec: contact — thermal contact.)
  bool has_thermal_contact() const;

  // Per-element active flag aligned with mesh.elements(): false for an id in
  // `deactivated_elements`, true otherwise. All-true when no element was removed, so
  // callers that gate their element loops on this mask are unchanged on a deck without
  // *MODEL CHANGE. (spec: model-change — deactivated elements carry no load.)
  std::vector<bool> element_active_mask() const;

  // Solver requested by the *STATIC step (SOLVER=), Auto when unspecified.
  RequestedSolver solver{RequestedSolver::Auto};

  // Number of eigenvalues requested by a *FREQUENCY step (first data-line field).
  // Inert unless procedure == Frequency. The eigensolution engine extracts the lowest
  // `num_eigenvalues` modes of K x = λ M x. (spec: modal-and-buckling — *FREQUENCY.)
  int num_eigenvalues{0};

  // Which complex-frequency eigenproblem a *COMPLEX FREQUENCY step requests. Only
  // Proportional is solved (the option-(B) proportional-damping reduction); Coriolis
  // (gyroscopic skew operator + rotor-speed body load) and Flutter (complex applied
  // force) are parsed then rejected — a faithful port needs input this deck does not
  // carry. (spec: modal-and-buckling-analysis — complex frequency.)
  enum class ComplexFrequencyType { Proportional, Coriolis, Flutter };

  // Number of complex modes requested by a *COMPLEX FREQUENCY step (data-line field).
  // Inert unless procedure == ComplexFrequency. The complex eigensolution engine reduces
  // onto the preceding *FREQUENCY basis and returns this many damped complex modes.
  int num_complex_modes{0};
  ComplexFrequencyType complex_freq_type{ComplexFrequencyType::Proportional};

  // Number of buckling modes requested by a *BUCKLE step (first data-line field).
  // Inert unless procedure == Buckling. The two-step prestress driver extracts the
  // lowest `num_buckling_modes` positive load factors of (K + λ K_geo) φ = 0. (spec:
  // input-deck-parsing — *BUCKLE; modal-and-buckling — *BUCKLE.)
  int num_buckling_modes{0};

  // Dynamics-step controls (spec: dynamic-analysis). Inert unless procedure is
  // ModalDynamic or SteadyStateDynamics. These carry the modal-superposition step
  // parameters parsed from *MODAL DYNAMIC / *STEADY STATE DYNAMICS and the damping model
  // (*DAMPING Rayleigh + *MODAL DAMPING ratios). The eigenbasis itself is extracted from
  // this same model (the mass/stiffness are procedure-independent), so a single-deck
  // modal-dynamic run reuses the *FREQUENCY extraction. `num_eigenvalues` selects how
  // many modes the superposition uses.
  Rayleigh rayleigh{};                 // *DAMPING, ALPHA=, BETA=
  std::vector<Real> modal_damping;     // *MODAL DAMPING ratios, indexed by mode
  Real dynamic_dt{0.0};                // *MODAL DYNAMIC time increment (data line 1)
  Real dynamic_t_end{0.0};             // *MODAL DYNAMIC step time (data line 2)
  Real steady_f_lo{0.0};               // *STEADY STATE DYNAMICS lower frequency
  Real steady_f_hi{0.0};               // *STEADY STATE DYNAMICS upper frequency
  int steady_num_points{0};            // *STEADY STATE DYNAMICS number of sweep points
  BaseMotion base_motion{};            // *BASE MOTION support excitation

  // *DYNAMIC direct-integration controls (spec: dynamic-analysis — direct dynamics).
  // Inert unless procedure == Dynamic. `dynamic_dt` / `dynamic_t_end` are shared with the
  // modal path (data line "dt, t_end"). `dynamic_alpha` is the HHT-α numerical-damping
  // knob (α ∈ [-1/3, 0]; 0 = energy-conserving trapezoidal Newmark), from the ALPHA=
  // parameter on *DYNAMIC. `dynamic_nonlinear` routes the step through the per-step Newton
  // driver (NLGEOM / a nonlinear material / contact), including the mass term in the
  // effective tangent; a linear model uses the factored effective-stiffness path.
  Real dynamic_alpha{0.0};             // *DYNAMIC, ALPHA= (HHT numerical damping)
  bool dynamic_nonlinear{false};       // *DYNAMIC, NLGEOM / nonlinear direct dynamics

  // *SUBSTRUCTURE GENERATE controls (spec: substructure-generation). Inert unless
  // procedure == Substructure. `retained_dofs` is the ordered set of retained (master)
  // nodal DOFs declared by *RETAINED NODAL DOFS — one entry per (node, comp) with comp
  // in 1..3. The exported superelement's DOFs are exactly these, in declaration order
  // (SORTED=NO) — the eliminated interior DOFs are Schur-condensed. The output flags
  // mirror *SUBSTRUCTURE MATRIX OUTPUT (STIFFNESS=/MASS=): a reduced mass is formed
  // (Craig-Bampton) only when `substructure_mass` is requested. `substructure_modes` is
  // the number of fixed-interface normal modes retained (0 -> pure Guyan static
  // reduction); it comes from a *FREQUENCY card inside the substructure step, or 0.
  struct RetainedDof {
    Index node_id{};
    int comp{};  // 1..3 (x,y,z)
  };
  std::vector<RetainedDof> retained_dofs;
  bool substructure_stiffness{true};   // *SUBSTRUCTURE MATRIX OUTPUT, STIFFNESS=YES
  bool substructure_mass{false};       // *SUBSTRUCTURE MATRIX OUTPUT, MASS=YES
  int substructure_modes{0};           // fixed-interface modes (Craig-Bampton); 0 = Guyan

  // Nonlinear-solution controls, parsed from *CONTROLS / *STATIC / *TIME POINTS.
  // Unused by the default linear path; consumed by solve_nonlinear_static.
  NonlinearControls controls{};
  Incrementation increment{};
  TimePoints time_points{};

  // Elastic properties per element (aligned with mesh.elements()), resolved from
  // the solid sections. Throws std::runtime_error on a missing elset/material or an
  // element left without a section.
  std::vector<ElasticIso> element_elastic() const;

  // Plastic properties per element (aligned with mesh.elements()), resolved from the
  // solid sections' materials. Elements whose material has no *PLASTIC data get an
  // empty optional (elastic-only). Last-writer wins per element, like element_elastic.
  std::vector<std::optional<Plastic>> element_plastic() const;

  // True when any element's material carries *PLASTIC data. A plastic model routes
  // the analysis to solve_nonlinear_static (load applied incrementally by the
  // driver); a purely elastic model keeps the linear path unchanged.
  bool has_plasticity() const;

  // Hyperelastic (*HYPERELASTIC) properties per element, aligned with
  // mesh.elements(); empty optional where the material is not hyperelastic.
  std::vector<std::optional<Hyperelastic>> element_hyperelastic() const;

  // User-material (*USER MATERIAL) properties per element, aligned with
  // mesh.elements(); empty optional where the material has no *USER MATERIAL.
  std::vector<std::optional<UserMaterial>> element_user_material() const;

  // True when any material carries a nonlinear constitutive law (*PLASTIC,
  // *HYPERELASTIC, or *USER MATERIAL). Routes the analysis to the incremental
  // nonlinear driver; a purely linear-elastic model keeps the linear path.
  bool has_nonlinear_material() const;

  // Mass density per element (aligned with mesh.elements()), resolved from the
  // solid sections' materials. Elements without *DENSITY get 0. Used by body loads.
  std::vector<Real> element_density() const;

  // Isotropic thermal conductivity per element (aligned with mesh.elements()),
  // resolved from the solid sections' materials. Throws on a missing elset/material
  // or an element left without a section (like element_elastic). An element whose
  // material carries no *CONDUCTIVITY throws — a heat-transfer solve needs it.
  // Uses the FIRST (constant) table value; the temperature-dependent evaluation is
  // element_conductivity_at(). (spec: heat-transfer-analysis — *CONDUCTIVITY.)
  std::vector<Real> element_conductivity() const;

  // True when any consumed thermal property (*CONDUCTIVITY / *SPECIFIC HEAT /
  // *EXPANSION) is a multi-row temperature table, so the thermal solve must
  // re-evaluate the property at the current temperature (Picard iteration). When
  // false every table is constant and the solve is byte-for-byte the pre-table path.
  bool has_temp_dependent_thermal() const;

  // Per-element conductivity evaluated at each element's MEAN nodal temperature from
  // `node_temp` (aligned with mesh node indices), for temperature-dependent
  // *CONDUCTIVITY k(T). Elements with a constant (single-row) table return that
  // constant regardless of `node_temp`, so a constant deck is unchanged.
  std::vector<Real> element_conductivity_at(const std::vector<Real>& node_temp) const;

  // Per-element volumetric heat capacity rho*c(T) evaluated at each element's mean
  // nodal temperature (temperature-dependent *SPECIFIC HEAT); constant tables are
  // unchanged. Aligned with mesh.elements().
  std::vector<Real> element_heat_capacity_at(const std::vector<Real>& node_temp) const;

  // Thermal-expansion data per element (aligned with mesh.elements()), resolved from
  // the solid sections' materials. Elements whose material carries no *EXPANSION get
  // an empty optional (no thermal strain). Last-writer wins per element, like
  // element_elastic. (spec: heat-transfer — thermal expansion coupling.)
  std::vector<std::optional<Expansion>> element_expansion() const;

  // True when any element's material carries *EXPANSION AND an applied temperature
  // field is present. When true a mechanical solve adds the thermal-strain load and
  // subtracts eps_th in stress recovery; when false the mechanical path is unchanged.
  bool has_thermal_strain() const;

  // Volumetric heat capacity rho*c per element (aligned with mesh.elements()) for
  // the transient capacitance matrix. Elements without both *DENSITY and *SPECIFIC
  // HEAT get 0 (steady-state does not use it). (spec: heat-transfer-analysis.)
  std::vector<Real> element_heat_capacity() const;

  // Expand every constraint (the raw *EQUATION relations plus the *MPC / *RIGID BODY
  // / *COUPLING / *TIE cards) into a flat list of linear Equations, ready for the
  // dependent-DOF elimination at assembly. Node/element sets and surfaces are
  // resolved against the mesh. (spec: constraints — all constraint types reduce to
  // linear equations eliminated during assembly.)
  std::vector<Equation> expand_constraints() const;

  // Scale factor for a load/BC at step fraction `lambda` in [0,1] (physical step
  // time t = lambda * increment.total). With no amplitude the default linear ramp
  // returns `lambda`; with an amplitude the curve is sampled at physical time `t`.
  // An unknown amplitude name falls back to the default ramp.
  Real amplitude_factor(const std::string& name, Real lambda) const;
};

}  // namespace cxpp
