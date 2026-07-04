#pragma once
#include <cstddef>
#include <vector>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/types.hpp"

// Design-sensitivity driver (spec: design-optimization — sensitivity core).
//
// Computes the gradient dObjective/dx of a design response with respect to coordinate
// design variables by the ADJOINT method, reusing the primal linear-static solution.
// Reference math (reimplemented, NEVER copied): CalculiX src/sensitivitys.f and
// src/objective_shapeener_tot.f (the sensi_coor path).
//
// For the reduced free/free linear-static system A u = b (A = K_ff, b the reduced
// external load), with a coordinate design variable x:
//   - Compliance response  g = bᵀu = uᵀA u  (self-adjoint):
//         dg/dx = 2 uᵀ (db/dx) − uᵀ (dA/dx) u        — reuses u, no adjoint solve.
//   - Nodal-displacement response  g = cᵀu  (general linear, c a unit selector):
//         solve the adjoint  A λ = c   (SAME operator A — reuses the primal
//                                       factorization / assembly),
//         dg/dx = λᵀ (db/dx − (dA/dx) u).
// The operator derivatives dA/dx, db/dx are obtained SEMI-ANALYTICALLY: the design
// coordinate is perturbed by ±h, the reduced system re-assembled, and A/b centrally
// differenced — the standard shape-sensitivity technique (CalculiX likewise perturbs
// the geometry per design variable).
namespace cxpp::numerics {

// The gradient of one design response over all coordinate design variables. `dgdx[i]`
// is dObjective/dx for model.design_variables[i]; `response_name` echoes the response.
struct SensitivityResult {
  std::string response_name;
  std::vector<Real> dgdx;   // one entry per design variable (aligned with the model)
  Real objective{0.0};      // the response value at the current design
};

// Full result of a *SENSITIVITY step: one SensitivityResult per *DESIGN RESPONSE.
struct SensitivityReport {
  std::vector<SensitivityResult> responses;
};

// Run the sensitivity analysis for `model` (procedure == Sensitivity). Solves the
// primal once, then computes each response's gradient over the coordinate design
// variables. Throws std::runtime_error if the model declares no design variables or an
// unsupported response. `h` is the semi-analytic coordinate perturbation.
SensitivityReport solve_sensitivity(const Model& model, Real h = 1e-4);

}  // namespace cxpp::numerics
