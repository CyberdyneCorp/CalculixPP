#pragma once
#include <optional>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/results.hpp"
#include "calculixpp/core/types.hpp"
#include "calculixpp/numerics/linear_static.hpp"

// Newton-Raphson nonlinear-static driver (spec: nonlinear-solution-control).
//
// This is an ADDED path: the default *STATIC solve stays linear
// (solve_linear_static). The driver reproduces the linear solve exactly on a linear
// model — with u=0 the first residual equals the external load, the tangent is the
// linear K, so du solves K du = f_ext and the second iteration has zero residual.
//
// The tangent is assemble_linear_static's reduced free/free system (for linear
// elasticity K_t == K). The internal force is fem::internal_force. Both are wired so
// a nonlinear tangent/stress can drop in later without changing the loop.
namespace cxpp::numerics {

// Per-solve options that are not carried on the Model. `line_search` is OFF by
// default (spec: Optional line search); it scales the Newton update to reduce the
// residual on difficult increments.
struct NonlinearOptions {
  bool line_search{false};
  // Force a solver kind (overrides the model's SOLVER=); empty -> resolve from the
  // model / system size, exactly like solve_linear_static.
  std::optional<SolverKind> forced{};
};

// Diagnostics from a nonlinear solve (increments taken, total Newton iterations,
// whether the full load factor 1.0 was reached).
struct NonlinearReport {
  int increments{0};
  int iterations{0};
  int cutbacks{0};
  bool converged{false};
  Real final_load_factor{0.0};
};

// Solve a static step with the Newton-Raphson driver, ramping the load factor from
// 0 to 1 through the incrementation engine. Returns the same StaticFields as
// solve_linear_static (displacement, stress, strain, reaction). Optionally reports
// solver statistics via `report`.
StaticFields solve_nonlinear_static(const Model& model,
                                    const NonlinearOptions& opts = {},
                                    NonlinearReport* report = nullptr);

}  // namespace cxpp::numerics
