#pragma once
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

#include "calculixpp/core/amplitude.hpp"
#include "calculixpp/core/connector.hpp"
#include "calculixpp/core/constraint.hpp"
#include "calculixpp/core/constraint_cards.hpp"
#include "calculixpp/core/material.hpp"
#include "calculixpp/core/mesh.hpp"
#include "calculixpp/core/section.hpp"
#include "calculixpp/core/solution_control.hpp"

namespace cxpp {

// Single-point constraint: prescribed nodal DOF (spec: loads-and-boundary-conditions).
// `amplitude` (optional *AMPLITUDE name) scales the prescribed value over step time.
struct Spc {
  Index node_id{};
  int comp{};      // 1..3 (x,y,z)
  Real value{0.0};
  std::string amplitude;  // empty -> default linear ramp
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

// The assembled analysis model for the linear-static slice.
class Model {
 public:
  Mesh mesh;
  std::unordered_map<std::string, Material> materials;
  std::vector<SolidSection> sections;
  std::vector<Spc> spcs;
  std::vector<Cload> cloads;
  std::vector<Dload> dloads;
  std::vector<BodyLoad> body_loads;
  std::unordered_map<std::string, Amplitude> amplitudes;

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

  // Solver requested by the *STATIC step (SOLVER=), Auto when unspecified.
  RequestedSolver solver{RequestedSolver::Auto};

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
