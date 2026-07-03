#include "calculixpp/numerics/multistep.hpp"

#include <stdexcept>

#include "calculixpp/fem/assembly.hpp"
#include "calculixpp/fem/loads.hpp"
#include "calculixpp/fem/stress.hpp"
#include "calculixpp/numerics/nonlinear_static.hpp"

namespace cxpp::numerics {
namespace {

// True when ANY step of the deck carries a *CONTACT PAIR — the deck is a multi-step
// contact analysis (spec: model-change — contact-pair activation between steps). Such a
// deck cannot use the linear-superposition incremental path (contact is a nonlinear
// constraint that does not superpose); it is solved per step through the nonlinear driver
// with each step's ACTIVE contact-pair set (with_active_contact_pairs), so a pair added or
// removed between steps changes that step's response.
bool any_contact(const std::vector<Model>& steps) {
  for (const Model& m : steps)
    if (m.has_contact()) return true;
  return false;
}

// Reject a step that this linear multi-step slice cannot yet handle. Multi-step
// thermal/coupled and nonlinear-material combinations are out of scope here — they
// need the thermal/nonlinear state carry, deferred with the contact/nonlinear
// multi-step follow-up. A clear throw beats a silently-wrong solve.
void require_linear_static(const Model& m) {
  if (m.procedure != Procedure::Static)
    throw std::runtime_error(
        "multi-step analysis currently supports only mechanical *STATIC steps; a "
        "thermal/coupled step in a multi-step deck is not implemented yet");
  if (m.has_nonlinear_material())
    throw std::runtime_error(
        "multi-step analysis currently supports only linear-elastic steps; a "
        "nonlinear-material (*PLASTIC/*HYPERELASTIC/*USER MATERIAL) multi-step deck "
        "is not implemented yet");
}

// Solve a multi-step CONTACT deck (spec: model-change — activate/deactivate a contact pair
// between steps). Each step is solved INDEPENDENTLY at its full applied load through the
// nonlinear contact driver, with its ACTIVE contact-pair set resolved from the accumulated
// *MODEL CHANGE, TYPE=CONTACT PAIR records (with_active_contact_pairs). Contact is a
// nonlinear constraint that does not superpose, so the incremental linear accumulation used
// for pure-elastic multi-step is not applicable here; solving each step's full model is the
// honest per-step model-change semantics (CalculiX re-solves the changed model each step).
// A step with no active contact pair falls back to the linear/nonlinear-material driver, so
// a pair inactive in one step and active in the next produces a different response — the
// required validation.
std::vector<StaticFields> solve_multistep_contact(const std::vector<Model>& steps,
                                                  std::optional<SolverKind> forced,
                                                  MultiStepReport* report) {
  std::vector<StaticFields> out;
  out.reserve(steps.size());
  for (std::size_t k = 0; k < steps.size(); ++k) {
    const Model step = steps[k].with_active_contact_pairs();
    if (step.procedure != Procedure::Static)
      throw std::runtime_error(
          "multi-step contact analysis supports only mechanical *STATIC steps");
    StaticFields f;
    if (step.has_contact() || step.has_nonlinear_material()) {
      NonlinearOptions opts;
      opts.forced = forced;
      f = solve_nonlinear_static(step, opts);
    } else {
      f = solve_linear_static(step, forced);
    }
    out.push_back(std::move(f));
    if (report) report->steps.push_back({static_cast<int>(k + 1), 0, step.procedure});
  }
  return out;
}

// Gather the free-DOF components of a full nodal displacement field into the reduced
// numbering of `sys` (u_prev_f[eq] = u_prev at that DOF). Constrained DOFs (SPC or MPC
// slave) have dof_eq < 0 and are skipped — their contribution is carried by the
// prescribed increment / the constraint transform, exactly as in a single-step solve.
std::vector<Real> gather_free(const fem::LinearSystem& sys,
                              const std::vector<Vec3>& u_prev) {
  std::vector<Real> uf(static_cast<std::size_t>(sys.n_free), 0.0);
  const std::size_t n_nodes = u_prev.size();
  for (std::size_t ni = 0; ni < n_nodes; ++ni)
    for (int c = 0; c < kDofsPerNode; ++c) {
      const std::size_t g = ni * kDofsPerNode + static_cast<std::size_t>(c);
      const Index eq = sys.dof_eq[g];
      if (eq >= 0)
        uf[static_cast<std::size_t>(eq)] = u_prev[ni][static_cast<std::size_t>(c)];
    }
  return uf;
}

// Reduced free/free matrix-vector product K_ff * x from the COO triplets.
std::vector<Real> spmv(const fem::LinearSystem& sys, const std::vector<Real>& x) {
  std::vector<Real> y(static_cast<std::size_t>(sys.n_free), 0.0);
  for (std::size_t k = 0; k < sys.vals.size(); ++k)
    y[static_cast<std::size_t>(sys.rows[k])] +=
        sys.vals[k] * x[static_cast<std::size_t>(sys.cols[k])];
  return y;
}

// Expand a reduced free-DOF vector to full nodal displacement through the step's
// constraint transform (free DOFs take uf[eq], SPC DOFs their prescribed value, MPC
// slaves reconstructed from masters) — the same expansion solve_linear_static uses.
std::vector<Vec3> expand(const fem::LinearSystem& sys, const std::vector<Real>& uf,
                         std::size_t n_nodes) {
  std::vector<Vec3> u(n_nodes, Vec3{0, 0, 0});
  for (std::size_t ni = 0; ni < n_nodes; ++ni)
    for (int c = 0; c < kDofsPerNode; ++c) {
      const std::size_t g = ni * kDofsPerNode + static_cast<std::size_t>(c);
      u[ni][static_cast<std::size_t>(c)] =
          sys.transform.displacement(g, uf, sys.prescribed);
    }
  return u;
}

// Solve one step for the displacement INCREMENT du (full nodal) from the carried total
// u_prev. K_ff du_f = sys.rhs - K_ff u_prev_f gives du_f such that u_prev+du is the
// step's total solution, while du is what a newly-active element sees (it enters K_ff
// this step, so it accrues strain only from du). Prescribed DOFs move by their increment
// (uc(k) - u_prev), reconstructed by expanding with prescribed - u_prev.
std::vector<Vec3> solve_step_increment(const fem::LinearSystem& sys,
                                       const std::vector<Vec3>& u_prev, SolverKind kind) {
  const std::size_t n_nodes = u_prev.size();
  const std::vector<Real> u_prev_f = gather_free(sys, u_prev);
  const std::vector<Real> k_uprev = spmv(sys, u_prev_f);
  fem::LinearSystem incr = sys;  // reuse the COO/transform; override the rhs
  for (std::size_t i = 0; i < incr.rhs.size(); ++i) incr.rhs[i] -= k_uprev[i];
  const std::vector<Real> duf = solve_reduced(incr, kind);

  std::vector<Real> prescribed_incr = sys.prescribed;
  for (std::size_t ni = 0; ni < n_nodes; ++ni)
    for (int c = 0; c < kDofsPerNode; ++c) {
      const std::size_t g = ni * kDofsPerNode + static_cast<std::size_t>(c);
      if (sys.transform.prescribed(g))
        prescribed_incr[g] -= u_prev[ni][static_cast<std::size_t>(c)];
    }
  fem::LinearSystem incr_sys = sys;
  incr_sys.prescribed = prescribed_incr;
  return expand(incr_sys, duf, n_nodes);
}

// Running accumulators for the incremental multi-step state (displacement total,
// stress/strain sums, internal-force sum). Each step advances them by its du and emits
// a StaticFields snapshot of the total state at that step's end.
struct StepAccum {
  std::vector<Vec3> u_total;
  std::vector<Voigt6> stress_sum;
  std::vector<Voigt6> strain_sum;
  std::vector<Real> f_int_total;

  explicit StepAccum(std::size_t n)
      : u_total(n, Vec3{0, 0, 0}),
        stress_sum(n, Voigt6{}),
        strain_sum(n, Voigt6{}),
        f_int_total(n * kDofsPerNode, 0.0) {}

  // Advance by the step increment du (stress/strain increment from recover_fields on
  // the active set, internal-force increment from internal_force) and return the total
  // reported field (reaction = accumulated f_int minus this step's external load).
  StaticFields advance(const Model& step, const std::vector<Vec3>& du) {
    const std::size_t n = u_total.size();
    StaticFields incr;
    incr.displacement = du;
    fem::recover_fields(step, incr);  // stress/strain from du on the step's active set
    const std::vector<Real> df_int = fem::internal_force(step, du);
    for (std::size_t i = 0; i < n; ++i) {
      for (int c = 0; c < kDofsPerNode; ++c)
        u_total[i][static_cast<std::size_t>(c)] += du[i][static_cast<std::size_t>(c)];
      for (int comp = 0; comp < 6; ++comp) {
        stress_sum[i][static_cast<std::size_t>(comp)] += incr.stress[i][static_cast<std::size_t>(comp)];
        strain_sum[i][static_cast<std::size_t>(comp)] += incr.strain[i][static_cast<std::size_t>(comp)];
      }
    }
    for (std::size_t i = 0; i < f_int_total.size(); ++i) f_int_total[i] += df_int[i];

    StaticFields sf;
    sf.displacement = u_total;
    sf.stress = stress_sum;
    sf.strain = strain_sum;
    sf.reaction.assign(n, Vec3{0, 0, 0});
    const std::vector<Real> f_ext = fem::external_load_vector(step);
    for (std::size_t i = 0; i < n; ++i)
      for (int c = 0; c < kDofsPerNode; ++c)
        sf.reaction[i][static_cast<std::size_t>(c)] =
            f_int_total[i * kDofsPerNode + static_cast<std::size_t>(c)] -
            f_ext[i * kDofsPerNode + static_cast<std::size_t>(c)];
    return sf;
  }
};

// Resolve any *BOUNDARY, FIXED SPC of `step` to the carried displacement u_prev at its
// DOF (hold the DOF at its deformed value across the boundary). Returns the step model
// unchanged when it has no FIXED SPC (the common case), so a deck without FIXED is
// untouched. In step 1 u_prev is zero, so a FIXED SPC prescribes 0 — the DOF's start
// state — which is the single-step interpretation.
Model resolve_fixed_bcs(const Model& step, const std::vector<Vec3>& u_prev) {
  bool any = false;
  for (const Spc& s : step.spcs)
    if (s.fixed) { any = true; break; }
  if (!any) return step;
  Model m = step;
  for (Spc& s : m.spcs) {
    if (!s.fixed) continue;
    if (s.comp < 1 || s.comp > kDofsPerNode) continue;
    const Index ni = m.mesh.node_index(s.node_id);
    s.value = u_prev[static_cast<std::size_t>(ni)][static_cast<std::size_t>(s.comp - 1)];
    s.amplitude.clear();  // a held (deformed) value is not ramped
  }
  return m;
}

}  // namespace

std::vector<StaticFields> solve_multistep_static_all(const std::vector<Model>& steps,
                                                     std::optional<SolverKind> forced,
                                                     MultiStepReport* report) {
  if (steps.empty()) return {};

  // A multi-step (or single-step) deck carrying any *CONTACT PAIR is a nonlinear-constraint
  // analysis: solve each step's model independently at full load with its active contact
  // pairs (contact does not superpose, so the linear incremental path is inapplicable).
  // (spec: model-change — activate/deactivate a contact pair between steps.)
  if (any_contact(steps)) return solve_multistep_contact(steps, forced, report);

  // Single-step fast path: EXACTLY the historical solve, so a one-*STEP deck is
  // byte-for-byte unchanged (the critical gate). No incremental machinery runs.
  if (steps.size() == 1) {
    require_linear_static(steps.front());
    StaticFields f = solve_linear_static(steps.front(), forced);
    if (report) report->steps.push_back({1, 0, steps.front().procedure});
    return {f};
  }

  const std::size_t n_nodes = steps.front().mesh.num_nodes();
  StepAccum acc(n_nodes);  // carried total displacement + stress/strain/f_int sums

  std::vector<StaticFields> out;
  out.reserve(steps.size());

  for (std::size_t k = 0; k < steps.size(); ++k) {
    require_linear_static(steps[k]);
    // Resolve *BOUNDARY, FIXED DOFs to their carried (deformed) values before assembly.
    const Model step = resolve_fixed_bcs(steps[k], acc.u_total);
    const fem::LinearSystem sys = fem::assemble_linear_static(step);
    const SolverKind kind =
        forced ? *forced : resolve_solver_kind(step.solver, sys.n_free);
    const std::vector<Vec3> du = solve_step_increment(sys, acc.u_total, kind);
    out.push_back(acc.advance(step, du));

    if (report)
      report->steps.push_back(
          {static_cast<int>(k + 1), sys.n_free, step.procedure});
  }

  return out;
}

StaticFields solve_multistep_static(const std::vector<Model>& steps,
                                    std::optional<SolverKind> forced,
                                    MultiStepReport* report) {
  std::vector<StaticFields> all = solve_multistep_static_all(steps, forced, report);
  if (all.empty()) return {};
  return std::move(all.back());
}

}  // namespace cxpp::numerics
