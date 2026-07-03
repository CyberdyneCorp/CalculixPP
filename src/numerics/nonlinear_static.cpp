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

// Reduce a full nodal DOF vector (a force/residual) to the free equations through the
// constraint transform (T^T f): a force on a free DOF lands on its equation; a force
// on an MPC slave DOF is distributed to its masters by the equation coefficients.
// Correct for both the residual (internal force) and the external load.
std::vector<Real> reduce_free(const fem::LinearSystem& sys,
                              const std::vector<Real>& full) {
  std::vector<Real> f(static_cast<std::size_t>(sys.n_free), 0.0);
  const fem::ConstraintTransform& tf = sys.transform;
  for (std::size_t g = 0; g < full.size(); ++g) {
    if (full[g] == 0.0) continue;
    for (const fem::DofTerm& t : tf.expansion[g].terms)
      f[static_cast<std::size_t>(t.eq)] += t.coeff * full[g];
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
  // Resolve the current (load-factor-scaled) SPC values per global DOF, then rebuild
  // free / SPC / MPC-slave motion through the constraint transform. Slaves are
  // reconstructed from their masters' current displacement, so an MPC that ties a
  // free DOF to a prescribed DOF tracks the amplitude ramp correctly.
  std::vector<Real> spc(sys.dof_eq.size(), 0.0);
  for (std::size_t g = 0; g < sys.dof_eq.size(); ++g)
    if (sys.transform.prescribed(g)) {
      const Real s = model.amplitude_factor(sys.prescribed_amp[g], lambda);
      spc[g] = s * sys.prescribed[g];
    }
  std::vector<Vec3> u(n_nodes, Vec3{0, 0, 0});
  for (std::size_t g = 0; g < sys.dof_eq.size(); ++g)
    u[g / kDofsPerNode][g % kDofsPerNode] = sys.transform.displacement(g, uf, spc);
  return u;
}

// Reduced external load at step fraction lambda: each load scaled by its amplitude
// (or the default linear ramp), reduced to the free DOFs.
std::vector<Real> external_free_at(const Model& model,
                                   const fem::LinearSystem& sys, Real lambda) {
  return reduce_free(sys, fem::external_load_vector(model, lambda));
}

// State carried across increments: the converged free-DOF solution, the DOF map /
// prescribed data (sys), the per-element material models + integration-point state,
// and the full-magnitude external-load norm for the relative force test.
struct NewtonState {
  fem::LinearSystem sys;          // DOF map + prescribed BCs (tangent reassembled)
  fem::MaterialPoints mp;         // per-element material models + Gauss-point state
  std::vector<Real> uf;           // current free-DOF solution
  Real f_ext_norm{0.0};           // ||f_ext,f|| at full magnitude
};

// Assemble the material-point tangent + internal force at free state uf and step
// fraction lambda. Returns the reduced tangent (with the DOF map from state.sys) and
// writes the reduced residual r = f_ext,f(lambda) - f_int,f into `r_out`.
fem::LinearSystem assemble_at(const Model& model, NewtonState& state,
                              const std::vector<Real>& uf, Real lambda,
                              std::vector<Real>& r_out,
                              Real* f_int_norm = nullptr) {
  const std::vector<Vec3> u =
      expand_displacement(model, state.sys, uf, lambda, model.mesh.num_nodes());
  std::vector<Real> f_int;
  fem::LinearSystem tangent =
      fem::assemble_material_tangent(model, u, state.mp, f_int);
  const std::vector<Real> f_int_free = reduce_free(tangent, f_int);
  const std::vector<Real> f_ext_free = external_free_at(model, tangent, lambda);
  r_out.assign(static_cast<std::size_t>(tangent.n_free), 0.0);
  for (std::size_t i = 0; i < r_out.size(); ++i)
    r_out[i] = f_ext_free[i] - f_int_free[i];
  if (f_int_norm) *f_int_norm = norm2(f_int_free);
  return tangent;
}

// Residual r = f_ext,f(lambda) - f_int,f on the free DOFs, at free state uf and step
// fraction lambda (recomputed from the material-point internal force). Used by the
// line search, which only needs the residual norm, not the tangent.
std::vector<Real> residual(const Model& model, NewtonState& state,
                           const std::vector<Real>& uf, Real lambda) {
  std::vector<Real> r;
  assemble_at(model, state, uf, lambda, r);
  return r;
}

// Line-search scale in (0,1]: pick alpha from a small set that most reduces the
// residual of the trial state uf + alpha*du (spec: Optional line search). OFF unless
// requested; returns 1.0 when disabled.
Real line_search_alpha(const Model& model, NewtonState& state,
                       const std::vector<Real>& uf, const std::vector<Real>& du,
                       Real lambda) {
  Real best_alpha = 1.0;
  Real best_norm = -1.0;
  for (Real alpha : {1.0, 0.5, 0.25, 0.125}) {
    std::vector<Real> trial(uf.size());
    for (std::size_t i = 0; i < uf.size(); ++i) trial[i] = uf[i] + alpha * du[i];
    const Real rn = norm2(residual(model, state, trial, lambda));
    if (best_norm < 0.0 || rn < best_norm) {
      best_norm = rn;
      best_alpha = alpha;
    }
  }
  return best_alpha;
}

// Solve the tangent system K_t du = r for the Newton correction.
std::vector<Real> solve_tangent(const Model& model, const NonlinearOptions& opts,
                                fem::LinearSystem tangent,
                                const std::vector<Real>& r) {
  tangent.rhs = r;  // reuse the COO tangent; swap in the residual rhs
  const SolverKind kind =
      opts.forced ? *opts.forced
                  : resolve_solver_kind(model.solver, tangent.n_free);
  return solve_reduced(tangent, kind);
}

// One Newton loop at load factor `lambda`, starting from state.uf. Returns the
// iteration count on convergence, or -1 if the iteration limit was hit (cutback).
// Each iteration reassembles the material-point tangent + internal-force residual at
// the current displacement, so nonlinear constitutive laws plug in unchanged.
// Convergence requires both the relative force residual and (after the first update)
// the relative displacement correction to fall below their tolerances.
int newton_iterate(const Model& model, const NonlinearOptions& opts,
                   NewtonState& state, Real lambda) {
  const NonlinearControls& ctl = model.controls;
  // Reference force for the relative residual test, fixed over the increment. For a
  // load-controlled step it is the applied load at this factor; for a displacement-
  // controlled step (||f_ext|| == 0) it is the reaction magnitude of the FIRST
  // iteration (the internal force that balances the prescribed motion). Holding it
  // fixed keeps ||r||/denom a genuine convergence measure — recomputing it each
  // iteration would make it identically 1 when f_ext == 0 (ref: CalculiX controlss.f
  // largest-force reference). (spec: nonlinear-solution-control — residual controls.)
  Real denom = 0.0;
  for (int it = 0; it < ctl.max_iterations; ++it) {
    std::vector<Real> r;
    Real f_int_norm = 0.0;
    const fem::LinearSystem tangent =
        assemble_at(model, state, state.uf, lambda, r, &f_int_norm);
    if (it == 0)
      denom = std::max({lambda * state.f_ext_norm, f_int_norm, ctl.eps});
    const std::vector<Real> du = solve_tangent(model, opts, tangent, r);
    const Real alpha =
        opts.line_search
            ? line_search_alpha(model, state, state.uf, du, lambda)
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
    // Increment converged: commit the per-integration-point trial history (plastic
    // strain / eqplastic / back stress) so the next increment return-maps from it.
    state.mp.commit();
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
  state.sys = fem::assemble_linear_static(model);  // DOF map + prescribed data
  state.mp = fem::make_material_points(model);      // per-element material + state
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
  // Plastic models recover stress/strain/reaction from the committed material-point
  // state (radial-return stress); elastic models keep the linear recovery unchanged.
  if (model.has_nonlinear_material())
    fem::recover_fields(model, res, state.mp);
  else
    fem::recover_fields(model, res);
  return res;
}

}  // namespace cxpp::numerics
