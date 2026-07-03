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
// One-way (sequential) coupling: solve the steady thermal field, copy the resulting
// nodal temperatures into the model's applied_temperature, then solve the mechanical
// field so the thermal strain eps_th = alpha (T - Tref) induces the thermal stress.
// Returns both fields. Full monolithic 4-DOF/node coupling is deferred; this
// captures the dominant (thermal->mechanical) coupling for thermal-stress problems.
CoupledFields solve_coupled(const Model& model,
                            std::optional<SolverKind> forced = std::nullopt);

}  // namespace cxpp::numerics
