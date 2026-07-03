#pragma once
#include <optional>
#include <vector>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/results.hpp"
#include "calculixpp/numerics/linear_static.hpp"

// Multi-step linear-static driver (spec: multi-step analysis; enables the cross-step
// deferrals — element birth-death 5.1, OP=MOD/NEW load accumulation 2.4, and *CHANGE
// SOLID SECTION / *CHANGE MATERIAL rebind 2.3). Each *STEP...*END STEP block is one
// Model (from io::parse_inp_steps); this driver solves them in order, carrying the
// converged displacement/stress/strain forward and accumulating them INCREMENTALLY so
// that:
//   - a constant active set reproduces the single-step total-load solve exactly
//     (linear superposition), so a two-load-step deck sums to the one-step result;
//   - an element removed in a later step contributes nothing thereafter (its active
//     mask is false), and one re-added is STRAIN-FREE relative to the DEFORMED
//     configuration at reactivation (it only accumulates strain from the increments of
//     the steps where it is active);
//   - a *BOUNDARY, ..., FIXED DOF holds its currently-attained (deformed) value across
//     the boundary (its prescribed increment is zero).
namespace cxpp::numerics {

// Per-step diagnostics from a multi-step solve.
struct StepReport {
  int index{0};        // 1-based step number
  Index n_free{0};     // free-DOF count of that step's reduced system
  Procedure procedure{Procedure::Static};
};

struct MultiStepReport {
  std::vector<StepReport> steps;
};

// Solve a list of *STEP models in order and return the final-step total fields
// (displacement/stress/strain/reaction), accumulated incrementally across steps.
// A single-element list is solved directly with solve_linear_static, so a one-*STEP
// deck is byte-for-byte the existing single-step result (the critical gate). Every
// step must be a mechanical (Static) linear-elastic step; a thermal/coupled or
// nonlinear-material step in a multi-step list throws (those multi-step combinations
// are out of this slice). `forced` overrides the solver kind. Optionally reports
// per-step diagnostics via `report`.
StaticFields solve_multistep_static(const std::vector<Model>& steps,
                                    std::optional<SolverKind> forced = std::nullopt,
                                    MultiStepReport* report = nullptr);

// Per-step total fields (index i is the accumulated state at the END of step i+1).
// Same solve as solve_multistep_static but returns every step's field, for
// per-step validation against CalculiX's per-step *NODE PRINT / *EL PRINT output.
std::vector<StaticFields> solve_multistep_static_all(
    const std::vector<Model>& steps, std::optional<SolverKind> forced = std::nullopt,
    MultiStepReport* report = nullptr);

}  // namespace cxpp::numerics
