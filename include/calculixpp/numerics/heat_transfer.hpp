#pragma once
#include <optional>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/results.hpp"
#include "calculixpp/numerics/linear_static.hpp"

// Steady-state heat-transfer solve (spec: heat-transfer-analysis). Assembles the
// scalar conduction system Kt T = q (via fem::assemble_conduction), applies the
// prescribed temperatures + fluxes, and solves the SPD system through the same
// ComputeBackend / SciPP path as the mechanical solver (numerics::solve_reduced).
// Returns nodal temperatures (NT) and the heat-flux reactions (RFL). This is a
// PARALLEL path — the mechanical solve is unchanged.
namespace cxpp::numerics {

// Solve a *HEAT TRANSFER step end to end, dispatching on the model's procedure:
// STEADY STATE solves Kt T = q (with film + iterated surface radiation); a
// transient step time-integrates (C/dt + Kt) T_{n+1} = (C/dt) T_n + q by backward
// Euler over the step period, returning the final-time field. `forced` overrides
// the solver kind (Direct by default; the operators are SPD). (spec: heat-transfer.)
ThermalFields solve_heat_transfer(const Model& model,
                                  std::optional<SolverKind> forced = std::nullopt);

// Solve a *COUPLED TEMPERATURE-DISPLACEMENT step (spec: heat-transfer — coupled).
// Dispatches on model.coupled_scheme:
//   - Monolithic: assemble+solve the 4-DOF/node (u,v,w,T) system in one linear solve.
//     The thermal-strain coupling enters as the off-diagonal K_uT block; the
//     mechanical->thermal block (plastic-dissipation heat, when taylor_quinney>0 with
//     plasticity) is added by an outer fixed point so the joint state converges.
//   - Staggered: Gauss-Seidel — alternate a thermal and a mechanical solve with an
//     outer convergence check; converges in one pass for one-way coupling.
// Both schemes reduce to the sequential thermal-then-mechanical solve for a pure
// thermal-stress problem (no mechanical->thermal feedback), so a thermal-stress deck
// gives the identical result under either scheme. Returns both fields.
CoupledFields solve_coupled(const Model& model,
                            std::optional<SolverKind> forced = std::nullopt);

// Diagnostics from a coupled solve: the fields plus the number of outer coupling
// iterations taken (1 for a one-way problem; more when the two-way plastic-dissipation
// heat is iterated to convergence). Used by tests and the reporting path.
struct CoupledReport {
  CoupledFields fields;
  int outer_iterations{0};
};
CoupledReport solve_coupled_reported(const Model& model,
                                     std::optional<SolverKind> forced = std::nullopt);

}  // namespace cxpp::numerics
