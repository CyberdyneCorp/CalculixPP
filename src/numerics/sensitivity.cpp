#include "calculixpp/numerics/sensitivity.hpp"

#include <stdexcept>
#include <unordered_map>

#include "calculixpp/fem/assembly.hpp"
#include "calculixpp/numerics/linear_static.hpp"

namespace cxpp::numerics {
namespace {

// Densify the reduced free/free COO stiffness of `sys` into a full-storage
// operator, so A·v (used by uᵀ dA/dx u and the adjoint residual) is cheap. The
// validatable design models are small, so an n_free × n_free dense operator is fine.
std::vector<Real> densify(const fem::LinearSystem& sys) {
  const std::size_t n = static_cast<std::size_t>(sys.n_free);
  std::vector<Real> A(n * n, 0.0);
  for (std::size_t k = 0; k < sys.vals.size(); ++k)
    A[static_cast<std::size_t>(sys.rows[k]) * n +
      static_cast<std::size_t>(sys.cols[k])] += sys.vals[k];
  return A;
}

// y = A v for the dense reduced operator (n × n row-major).
std::vector<Real> matvec(const std::vector<Real>& A, const std::vector<Real>& v) {
  const std::size_t n = v.size();
  std::vector<Real> y(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    Real s = 0.0;
    for (std::size_t j = 0; j < n; ++j) s += A[i * n + j] * v[j];
    y[i] = s;
  }
  return y;
}

Real dot(const std::vector<Real>& a, const std::vector<Real>& b) {
  Real s = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
  return s;
}

// Assemble the reduced system for a model with one design coordinate perturbed by
// `delta` (0 => the base model). Returns the reduced operator A (dense) and load b on
// the SAME free-DOF numbering as the base system (the design perturbation is a small
// coordinate move that does not change the constraint topology, so the numbering is
// stable). `var` selects the (node, comp) coordinate.
struct ReducedSystem {
  std::vector<Real> A;  // n_free × n_free dense
  std::vector<Real> b;  // reduced load, size n_free
};

ReducedSystem assemble_perturbed(const Model& base, const DesignVariable& var,
                                 Real delta) {
  Model m = base;  // value copy — mutate its mesh coordinate only
  const Index ni = m.mesh.node_index(var.node_id);
  if (ni < 0) throw std::runtime_error("*DESIGN VARIABLES references unknown node");
  if (delta != 0.0) m.mesh.perturb_node_coord(ni, var.comp - 1, delta);
  const fem::LinearSystem sys = fem::assemble_linear_static(m);
  return {densify(sys), sys.rhs};
}

// The free-equation index of a nodal DOF (node id, comp 1..3), or -1 if the DOF is
// constrained. Uses the base system's dof_eq map.
Index free_equation(const Model& model, const fem::LinearSystem& sys, Index node_id,
                    int comp) {
  const Index ni = model.mesh.node_index(node_id);
  if (ni < 0) throw std::runtime_error("*DESIGN RESPONSE references unknown node");
  return sys.dof_eq[static_cast<std::size_t>(ni) * kDofsPerNode +
                    static_cast<std::size_t>(comp - 1)];
}

// Gradient of one response over every design variable, given the primal free solution
// `u`, the base reduced operator `A0`, and (for the displacement response) the adjoint
// multiplier `lambda`. For each variable we central-difference the reduced operator:
//   dA/dx = (A+ − A−)/2h,  db/dx = (b+ − b−)/2h,
// then apply the closed-form adjoint gradient.
std::vector<Real> response_gradient(const Model& model, const DesignResponse& resp,
                                    const fem::LinearSystem& base,
                                    const std::vector<Real>& u,
                                    const std::vector<Real>& b0,
                                    const std::vector<Real>& lambda, Real h) {
  const bool compliance = resp.kind == DesignResponse::Kind::Compliance;
  std::vector<Real> dgdx(model.design_variables.size(), 0.0);
  for (std::size_t v = 0; v < model.design_variables.size(); ++v) {
    const ReducedSystem sp = assemble_perturbed(model, model.design_variables[v], h);
    const ReducedSystem sm = assemble_perturbed(model, model.design_variables[v], -h);
    const std::size_t n = u.size();
    // dA/dx · u  and  db/dx
    std::vector<Real> dAu(n, 0.0), dbdx(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) dbdx[i] = (sp.b[i] - sm.b[i]) / (2.0 * h);
    // (dA/dx) u = ((A+ − A−)/2h) u
    const std::vector<Real> Apu = matvec(sp.A, u);
    const std::vector<Real> Amu = matvec(sm.A, u);
    for (std::size_t i = 0; i < n; ++i) dAu[i] = (Apu[i] - Amu[i]) / (2.0 * h);

    if (compliance) {
      // dg/dx = 2 uᵀ db/dx − uᵀ (dA/dx) u
      dgdx[v] = 2.0 * dot(u, dbdx) - dot(u, dAu);
    } else {
      // dg/dx = λᵀ (db/dx − (dA/dx) u)
      std::vector<Real> rhs(n, 0.0);
      for (std::size_t i = 0; i < n; ++i) rhs[i] = dbdx[i] - dAu[i];
      dgdx[v] = dot(lambda, rhs);
    }
  }
  (void)base;
  return dgdx;
}

}  // namespace

SensitivityReport solve_sensitivity(const Model& model, Real h) {
  if (model.design_variables.empty())
    throw std::runtime_error("*SENSITIVITY step has no *DESIGN VARIABLES");
  if (model.design_responses.empty())
    throw std::runtime_error("*SENSITIVITY step has no *DESIGN RESPONSE / *OBJECTIVE");

  // Primal: assemble and solve K u = f once. The reduced operator A0 is reused for
  // both the compliance objective (uᵀA0u) and the adjoint response solve (A0 λ = c).
  const fem::LinearSystem base = fem::assemble_linear_static(model);
  const std::vector<Real> u = solve_reduced(base, SolverKind::Direct);
  const std::vector<Real>& b0 = base.rhs;

  SensitivityReport report;
  for (const DesignResponse& resp : model.design_responses) {
    std::vector<Real> lambda;  // adjoint multiplier (displacement response only)
    Real objective = 0.0;
    if (resp.kind == DesignResponse::Kind::Compliance) {
      objective = dot(b0, u);  // g = fᵀu (external work / compliance)
    } else {
      const Index eq = free_equation(model, base, resp.node_id, resp.comp);
      if (eq < 0)
        throw std::runtime_error(
            "*DESIGN RESPONSE displacement DOF is constrained (zero sensitivity)");
      objective = u[static_cast<std::size_t>(eq)];
      // Adjoint A λ = c with c the unit selector of the response DOF. Reuses the SAME
      // assembled reduced operator as the primal (not a re-factor per design variable).
      std::vector<Real> c(u.size(), 0.0);
      c[static_cast<std::size_t>(eq)] = 1.0;
      fem::LinearSystem adj = base;  // same COO operator, different rhs
      adj.rhs = c;
      lambda = solve_reduced(adj, SolverKind::Direct);
    }
    SensitivityResult r;
    r.response_name = resp.name;
    r.objective = objective;
    r.dgdx = response_gradient(model, resp, base, u, b0, lambda, h);
    report.responses.push_back(std::move(r));
  }
  return report;
}

}  // namespace cxpp::numerics
