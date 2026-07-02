#include "calculixpp/numerics/nonlinear_static.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "calculixpp/fem/assembly.hpp"
#include "calculixpp/fem/loads.hpp"
#include "calculixpp/fem/stress.hpp"

namespace cxpp::numerics {
namespace {

// Euclidean norm of a reduced vector.
Real norm2(const std::vector<Real>& v) {
  Real s = 0.0;
  for (Real x : v) s += x * x;
  return std::sqrt(s);
}

// Reduce a full nodal DOF vector to the free equations using the DOF map.
std::vector<Real> reduce_free(const fem::LinearSystem& sys,
                              const std::vector<Real>& full) {
  std::vector<Real> f(static_cast<std::size_t>(sys.n_free), 0.0);
  for (std::size_t g = 0; g < sys.dof_eq.size(); ++g) {
    const Index eq = sys.dof_eq[g];
    if (eq >= 0) f[static_cast<std::size_t>(eq)] = full[g];
  }
  return f;
}

// Expand the free-DOF solution u_f and the prescribed BCs to a full nodal
// displacement field, ready for fem::internal_force. Prescribed DOFs are scaled by
// their *AMPLITUDE factor at step fraction lambda (default linear ramp when no
// amplitude is attached).
std::vector<Vec3> expand_displacement(const Model& model,
                                      const fem::LinearSystem& sys,
                                      const std::vector<Real>& uf, Real lambda,
                                      std::size_t n_nodes) {
  std::vector<Vec3> u(n_nodes, Vec3{0, 0, 0});
  for (std::size_t g = 0; g < sys.dof_eq.size(); ++g) {
    const std::size_t ni = g / kDofsPerNode;
    const std::size_t c = g % kDofsPerNode;
    const Index eq = sys.dof_eq[g];
    if (eq >= 0) {
      u[ni][c] = uf[static_cast<std::size_t>(eq)];
    } else {
      const Real s = model.amplitude_factor(sys.prescribed_amp[g], lambda);
      u[ni][c] = s * sys.prescribed[g];
    }
  }
  return u;
}

// Reduced external load at step fraction lambda: each load scaled by its amplitude
// (or the default linear ramp), reduced to the free DOFs.
std::vector<Real> external_free_at(const Model& model,
                                   const fem::LinearSystem& sys, Real lambda) {
  return reduce_free(sys, fem::external_load_vector(model, lambda));
}

// Residual r = f_ext,f(lambda) - f_int,f on the free DOFs, at free state uf and step
// fraction lambda. f_ext is amplitude-scaled per load; f_int is the global internal
// force (fem::internal_force), reduced.
std::vector<Real> residual(const Model& model, const fem::LinearSystem& sys,
                           const std::vector<Real>& uf, Real lambda) {
  const std::vector<Vec3> u =
      expand_displacement(model, sys, uf, lambda, model.mesh.num_nodes());
  const std::vector<Real> f_int_free =
      reduce_free(sys, fem::internal_force(model, u));
  const std::vector<Real> f_ext_free = external_free_at(model, sys, lambda);
  std::vector<Real> r(static_cast<std::size_t>(sys.n_free), 0.0);
  for (std::size_t i = 0; i < r.size(); ++i)
    r[i] = f_ext_free[i] - f_int_free[i];
  return r;
}

// Line-search scale in (0,1]: pick alpha from a small set that most reduces the
// residual of the trial state uf + alpha*du (spec: Optional line search). OFF unless
// requested; returns 1.0 when disabled.
Real line_search_alpha(const Model& model, const fem::LinearSystem& sys,
                       const std::vector<Real>& uf, const std::vector<Real>& du,
                       Real lambda) {
  Real best_alpha = 1.0;
  Real best_norm = -1.0;
  for (Real alpha : {1.0, 0.5, 0.25, 0.125}) {
    std::vector<Real> trial(uf.size());
    for (std::size_t i = 0; i < uf.size(); ++i) trial[i] = uf[i] + alpha * du[i];
    const Real rn = norm2(residual(model, sys, trial, lambda));
    if (best_norm < 0.0 || rn < best_norm) {
      best_norm = rn;
      best_alpha = alpha;
    }
  }
  return best_alpha;
}

// State carried across increments: the converged free-DOF solution and the assembled
// reduced system (tangent + external load), computed once for the linear path.
struct NewtonState {
  fem::LinearSystem sys;
  std::vector<Real> uf;   // current free-DOF solution
  Real f_ext_norm{0.0};   // ||f_ext,f|| at full magnitude, for the relative test
};

// Solve the tangent system K_t du = r for the Newton correction.
std::vector<Real> solve_tangent(const Model& model, const NonlinearOptions& opts,
                                const fem::LinearSystem& sys,
                                const std::vector<Real>& r) {
  fem::LinearSystem tangent = sys;  // reuse the COO tangent; swap in the residual rhs
  tangent.rhs = r;
  const SolverKind kind =
      opts.forced ? *opts.forced
                  : resolve_solver_kind(model.solver, tangent.n_free);
  return solve_reduced(tangent, kind);
}

// One Newton loop at load factor `lambda`, starting from state.uf. Returns the
// iteration count on convergence, or -1 if the iteration limit was hit (cutback).
// Convergence requires both the relative force residual and (after the first update)
// the relative displacement correction to fall below their tolerances.
int newton_iterate(const Model& model, const NonlinearOptions& opts,
                   NewtonState& state, Real lambda) {
  const NonlinearControls& ctl = model.controls;
  const Real denom = std::max(lambda * state.f_ext_norm, ctl.eps);
  for (int it = 0; it < ctl.max_iterations; ++it) {
    const std::vector<Real> r =
        residual(model, state.sys, state.uf, lambda);
    const std::vector<Real> du = solve_tangent(model, opts, state.sys, r);
    const Real alpha =
        opts.line_search
            ? line_search_alpha(model, state.sys, state.uf, du, lambda)
            : 1.0;
    Real du_norm = 0.0, u_norm = 0.0;
    for (std::size_t i = 0; i < state.uf.size(); ++i) {
      state.uf[i] += alpha * du[i];
      du_norm += (alpha * du[i]) * (alpha * du[i]);
      u_norm += state.uf[i] * state.uf[i];
    }
    // Converged when the force residual AND the displacement correction that was
    // just applied are both within tolerance.
    const Real disp_ratio =
        std::sqrt(du_norm) / std::max(std::sqrt(u_norm), ctl.eps);
    if (norm2(r) / denom < ctl.force_tol && disp_ratio < ctl.disp_tol)
      return it + 1;
  }
  return -1;  // did not converge -> cutback
}

// Clamp an increment so it lands on the next *TIME POINTS mark and never overshoots
// the end of the step (lambda == 1). Returns the (possibly reduced) increment.
Real clamp_increment(const Model& model, Real lambda, Real dl) {
  Real next = std::min<Real>(1.0, lambda + dl);
  for (Real t : model.time_points.times) {
    const Real f = t / std::max(model.increment.total, model.controls.eps);
    if (f > lambda + model.controls.eps && f < next) next = f;
  }
  return next - lambda;
}

// Drive the load factor 0 -> 1 with automatic incrementation and cutback. Updates
// state.uf in place; fills the report. Returns true if lambda reached 1.
bool run_increments(const Model& model, const NonlinearOptions& opts,
                    NewtonState& state, NonlinearReport& rep) {
  const Incrementation& inc = model.increment;
  const Real total = std::max(inc.total, model.controls.eps);
  Real lambda = 0.0;
  Real dl = std::clamp(inc.initial / total, inc.min / total, inc.max / total);
  while (lambda < 1.0 - model.controls.eps &&
         rep.increments < inc.max_increments) {
    const Real step = clamp_increment(model, lambda, dl);
    const std::vector<Real> saved = state.uf;
    const int iters = newton_iterate(model, opts, state, lambda + step);
    if (iters < 0) {  // cutback
      if (inc.direct || step * total <= inc.min + model.controls.eps) return false;
      state.uf = saved;
      dl = std::max(inc.min / total, dl * inc.cutback);
      ++rep.cutbacks;
      continue;
    }
    lambda += step;
    ++rep.increments;
    rep.iterations += iters;
    if (!inc.direct && iters <= inc.grow_below)
      dl = std::min(inc.max / total, dl * inc.grow);
  }
  rep.final_load_factor = lambda;
  return lambda >= 1.0 - model.controls.eps;
}

}  // namespace

StaticFields solve_nonlinear_static(const Model& model,
                                    const NonlinearOptions& opts,
                                    NonlinearReport* report) {
  NewtonState state;
  state.sys = fem::assemble_linear_static(model);
  state.f_ext_norm =
      norm2(reduce_free(state.sys, fem::external_load_vector(model)));
  state.uf.assign(static_cast<std::size_t>(state.sys.n_free), 0.0);

  NonlinearReport rep;
  rep.converged = run_increments(model, opts, state, rep);
  if (report) *report = rep;

  StaticFields res;
  const std::size_t n_nodes = model.mesh.num_nodes();
  res.displacement.assign(n_nodes, Vec3{0, 0, 0});
  const std::vector<Vec3> u = expand_displacement(
      model, state.sys, state.uf, rep.final_load_factor, n_nodes);
  res.displacement = u;
  fem::recover_fields(model, res);  // stress, strain, reaction
  return res;
}

}  // namespace cxpp::numerics
