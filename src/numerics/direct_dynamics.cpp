#include "calculixpp/numerics/direct_dynamics.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>

#include "calculixpp/fem/contact.hpp"
#include "calculixpp/fem/loads.hpp"
#include "calculixpp/fem/stress.hpp"

// Direct-integration dynamics (spec: dynamic-analysis; tasks 3.1/3.2/3.3). HHT-α
// implicit integration of M a + C v + K u = f(t) in physical coordinates. The linear
// path forms the effective dynamic stiffness once and re-solves each step through
// solve_reduced; the nonlinear path (task 3.2) re-assembles the material-point tangent
// each Newton iteration and adds the inertial/damping terms to it. Energy is tracked
// each step (task 3.3): an undamped α=0 run conserves total energy over the history.
namespace cxpp::numerics {
namespace {

using fem::LinearSystem;

// Dense symmetric n×n buffer from the free/free COO triplets of a LinearSystem
// (duplicates summed, both triangles filled). Same congruence-scatter symmetrization as
// the eigensolution path — the assembler stores each entry once (upper OR lower).
std::vector<Real> dense_symmetric(const LinearSystem& sys) {
  const std::size_t n = static_cast<std::size_t>(sys.n_free);
  std::vector<Real> A(n * n, 0.0);
  for (std::size_t k = 0; k < sys.vals.size(); ++k) {
    const std::size_t r = static_cast<std::size_t>(sys.rows[k]);
    const std::size_t c = static_cast<std::size_t>(sys.cols[k]);
    A[r * n + c] += sys.vals[k];
  }
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = i + 1; j < n; ++j) {
      const Real aij = A[i * n + j], aji = A[j * n + i];
      const Real val = (aij != 0.0 && aji != 0.0) ? aij : aij + aji;
      A[i * n + j] = val;
      A[j * n + i] = val;
    }
  return A;
}

// y = A x for a dense n×n symmetric buffer.
std::vector<Real> matvec(const std::vector<Real>& A, const std::vector<Real>& x,
                         std::size_t n) {
  std::vector<Real> y(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    Real s = 0.0;
    const Real* row = &A[i * n];
    for (std::size_t j = 0; j < n; ++j) s += row[j] * x[j];
    y[i] = s;
  }
  return y;
}

Real dot(const std::vector<Real>& a, const std::vector<Real>& b) {
  Real s = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
  return s;
}

// Reduce a full nodal DOF vector (a force) to the free equations through the constraint
// transform (Tᵀ f). Mirrors the nonlinear driver's reduce_free.
std::vector<Real> reduce_free(const LinearSystem& sys, const std::vector<Real>& full) {
  std::vector<Real> f(static_cast<std::size_t>(sys.n_free), 0.0);
  const fem::ConstraintTransform& tf = sys.transform;
  for (std::size_t g = 0; g < full.size(); ++g) {
    if (full[g] == 0.0) continue;
    for (const fem::DofTerm& t : tf.expansion[g].terms)
      f[static_cast<std::size_t>(t.eq)] += t.coeff * full[g];
  }
  return f;
}

// Expand a free-DOF vector to a full nodal field through the constraint transform (SPC
// value 0 — a dynamic response about the rest state has no prescribed motion here).
std::vector<Vec3> expand_free(const LinearSystem& sys, const std::vector<Real>& q,
                              std::size_t n_nodes) {
  const std::vector<Real> zero_prescribed(sys.prescribed.size(), 0.0);
  std::vector<Vec3> u(n_nodes, Vec3{0, 0, 0});
  for (std::size_t ni = 0; ni < n_nodes; ++ni)
    for (int c = 0; c < kDofsPerNode; ++c) {
      const std::size_t g = ni * kDofsPerNode + static_cast<std::size_t>(c);
      u[ni][static_cast<std::size_t>(c)] =
          sys.transform.displacement(g, q, zero_prescribed);
    }
  return u;
}

// The Newmark/HHT integration state on the free DOFs: displacement u, velocity v,
// acceleration a. Advanced one step by the HHT solver.
struct MotionState {
  std::vector<Real> u, v, a;
};

// HHT-α time-marching coefficients derived from (Δt, β, γ, α), grouped so the linear and
// nonlinear paths share them. Solving for the acceleration a_{n+1} directly (Newmark
// correctors substituted into the HHT equilibrium), the effective dynamic operator is
//     A_eff = M + (1+α)γΔt C + (1+α)βΔt² K = c_m M + c_c C + k_scale K,
// with c_m = 1, c_c = (1+α)γΔt, k_scale = (1+α)βΔt². The right-hand side is then
//     (1+α)f_{n+1} - αf_n - (1+α)(K ũ + C ṽ) + α(K uₙ + C vₙ),
// undivided (both sides in the same scaling). Predictor (a_{n+1}=0):
//     ũ = u + Δt v + Δt²(½-β) a,   ṽ = v + Δt(1-γ) a.
struct HhtCoeffs {
  Real dt, alpha, beta, gamma;
  Real c_m, c_c, k_scale;  // A_eff = c_m M + c_c C + k_scale K

  HhtCoeffs(Real dt_, const HhtParams& p)
      : dt(dt_), alpha(p.alpha), beta(p.beta), gamma(p.gamma) {
    c_m = 1.0;
    c_c = (1.0 + alpha) * gamma * dt;
    k_scale = (1.0 + alpha) * beta * dt * dt;
  }
};

// Scatter c_m*M + c_c*C + k_scale*K into a fresh COO effective-stiffness system, reusing
// `base` (K) for its DOF map / transform. C = damp.alpha*M + damp.beta*K. Because M, K
// share the free-DOF numbering, we accumulate a dense effective buffer then emit its
// triplets — small validation meshes make the dense path exact and simple.
LinearSystem effective_system(const LinearSystem& K, const std::vector<Real>& Md,
                              const std::vector<Real>& Kd, const HhtCoeffs& hc,
                              const Damping& damp, std::size_t n) {
  // C = alpha M + beta K  ->  c_c*C + c_m*M = (c_c*alpha + c_m) M + c_c*beta K.
  const Real m_coeff = hc.c_m + hc.c_c * damp.alpha;
  const Real k_coeff = hc.k_scale + hc.c_c * damp.beta;
  LinearSystem eff = K;
  eff.rows.clear();
  eff.cols.clear();
  eff.vals.clear();
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j) {
      const Real e = m_coeff * Md[i * n + j] + k_coeff * Kd[i * n + j];
      if (e == 0.0) continue;
      eff.rows.push_back(static_cast<Index>(i));
      eff.cols.push_back(static_cast<Index>(j));
      eff.vals.push_back(e);
    }
  return eff;
}

// Kinetic + strain energy of the current state (task 3.3). KE = ½ vᵀ M v; strain
// energy = ½ uᵀ K u (linear elastic). Total drift over an undamped run measures energy
// conservation.
void fill_energy(DirectTimePoint& tp, const std::vector<Real>& Md,
                 const std::vector<Real>& Kd, const MotionState& s, std::size_t n) {
  tp.kinetic_energy = 0.5 * dot(s.v, matvec(Md, s.v, n));
  tp.strain_energy = 0.5 * dot(s.u, matvec(Kd, s.u, n));
  tp.total_energy = tp.kinetic_energy + tp.strain_energy;
}

DirectTimePoint make_frame(const LinearSystem& K, const std::vector<Real>& Md,
                           const std::vector<Real>& Kd, const MotionState& s, Real t,
                           std::size_t n, std::size_t n_nodes) {
  DirectTimePoint tp;
  tp.time = t;
  tp.displacement = expand_free(K, s.u, n_nodes);
  tp.velocity = expand_free(K, s.v, n_nodes);
  tp.acceleration = expand_free(K, s.a, n_nodes);
  fill_energy(tp, Md, Kd, s, n);
  return tp;
}

// External free-DOF load at time t: spatial pattern (full-magnitude external load) scaled
// by amplitude(t) (default constant step load), reduced to free DOFs.
std::vector<Real> load_at(const std::vector<Real>& f_free,
                          const std::function<Real(Real)>& amp, Real t) {
  const Real s = amp ? amp(t) : 1.0;
  std::vector<Real> f(f_free.size());
  for (std::size_t i = 0; i < f.size(); ++i) f[i] = s * f_free[i];
  return f;
}

// Advance one HHT step of the LINEAR system. Solves K_eff Δ = R for the new acceleration
// increment via the effective-stiffness residual, then applies the Newmark corrector.
// The HHT residual solved for a_{n+1} directly (u,v expressed via a_{n+1}) is
//     K_eff a_{n+1} = (1+α) f_{n+1} - α f_n
//                     - (1+α) K ũ - c_c C-terms... — assembled below in free DOFs.
void step_linear(LinearSystem eff, const std::vector<Real>& Kd,
                 const std::vector<Real>& Cd, const HhtCoeffs& hc, MotionState& s,
                 const std::vector<Real>& f_new, const std::vector<Real>& f_old,
                 SolverKind kind, std::size_t n) {
  const Real dt = hc.dt, beta = hc.beta, gamma = hc.gamma, alpha = hc.alpha;
  // Predictors (state at a_{n+1}=0): u_p = u + dt v + dt²(½-β) a; v_p = v + dt(1-γ) a.
  std::vector<Real> u_p(n), v_p(n);
  for (std::size_t i = 0; i < n; ++i) {
    u_p[i] = s.u[i] + dt * s.v[i] + dt * dt * (0.5 - beta) * s.a[i];
    v_p[i] = s.v[i] + dt * (1.0 - gamma) * s.a[i];
  }
  // Effective load: (1+α) f_{n+1} - α f_n - (1+α)(K u_p + C v_p) + α (K u_n + C v_n).
  const std::vector<Real> Ku_p = matvec(Kd, u_p, n), Cv_p = matvec(Cd, v_p, n);
  const std::vector<Real> Ku_n = matvec(Kd, s.u, n), Cv_n = matvec(Cd, s.v, n);
  std::vector<Real> rhs(n);
  for (std::size_t i = 0; i < n; ++i)
    rhs[i] = (1.0 + alpha) * f_new[i] - alpha * f_old[i] -
             (1.0 + alpha) * (Ku_p[i] + Cv_p[i]) + alpha * (Ku_n[i] + Cv_n[i]);
  eff.rhs = rhs;
  const std::vector<Real> a_new = solve_reduced(eff, kind);
  // Correctors: u_{n+1} = u_p + β dt² a; v_{n+1} = v_p + γ dt a.
  for (std::size_t i = 0; i < n; ++i) {
    s.u[i] = u_p[i] + beta * dt * dt * a_new[i];
    s.v[i] = v_p[i] + gamma * dt * a_new[i];
    s.a[i] = a_new[i];
  }
}

// Newmark correctors that map a candidate acceleration a_{n+1} to u_{n+1}, v_{n+1} from
// the predictors u_p, v_p: u = u_p + β Δt² a; v = v_p + γ Δt a.
void newmark_correct(MotionState& s, const std::vector<Real>& u_p,
                     const std::vector<Real>& v_p, const std::vector<Real>& a,
                     const HhtCoeffs& hc, std::size_t n) {
  for (std::size_t i = 0; i < n; ++i) {
    s.u[i] = u_p[i] + hc.beta * hc.dt * hc.dt * a[i];
    s.v[i] = v_p[i] + hc.gamma * hc.dt * a[i];
    s.a[i] = a[i];
  }
}

// ---------------------------------------------------------------------------
// Nonlinear direct dynamics (task 3.2): the effective tangent per Newton iteration
// includes the mass/damping terms, and the residual carries inertia + internal force
// (+ contact). Driven through the material-point tangent so plasticity/hyperelasticity
// and contact ride along with inertia.
// ---------------------------------------------------------------------------

// Euclidean norm of a reduced vector.
Real norm2(const std::vector<Real>& v) { return std::sqrt(dot(v, v)); }

// State carried across the nonlinear Newton dynamics run.
struct NonlinDyn {
  fem::MaterialPoints mp;
  std::vector<fem::ResolvedContactPair> contact;
};

// Assemble the material-point tangent K_tan + internal force f_int at free displacement
// `uf`, returning the reduced tangent (COO) and writing the reduced internal force into
// `fint_free`. Mirrors the nonlinear-static assemble_at without the external load.
LinearSystem assemble_tangent(const Model& model, NonlinDyn& st,
                              const LinearSystem& K, const std::vector<Real>& uf,
                              std::size_t n_nodes, std::vector<Real>& fint_free) {
  const std::vector<Real> zero_prescribed(K.prescribed.size(), 0.0);
  std::vector<Vec3> u(n_nodes, Vec3{0, 0, 0});
  for (std::size_t g = 0; g < K.dof_eq.size(); ++g)
    u[g / kDofsPerNode][g % kDofsPerNode] =
        K.transform.displacement(g, uf, zero_prescribed);
  std::vector<Real> f_int;
  LinearSystem tangent = fem::assemble_material_tangent(model, u, st.mp, f_int);
  if (!st.contact.empty())
    fem::add_contact(model, st.contact, u, tangent, f_int);
  fint_free = reduce_free(tangent, f_int);
  return tangent;
}

// Fold the material tangent into the effective dynamic tangent, in place:
//     K_eff = c_m M + c_c C + k_scale K_tan,   C = α M + β K_lin (dense Cd).
// The material tangent COO is densified, scaled by k_scale, has the inertia/damping
// dense buffers added, and re-emitted as summed symmetric triplets.
void fold_effective_tangent(LinearSystem& tangent, const std::vector<Real>& Md,
                            const std::vector<Real>& Cd, const HhtCoeffs& hc,
                            std::size_t n) {
  std::vector<Real> T = dense_symmetric(tangent);
  for (std::size_t k = 0; k < n * n; ++k)
    T[k] = hc.k_scale * T[k] + hc.c_m * Md[k] + hc.c_c * Cd[k];
  tangent.rows.clear();
  tangent.cols.clear();
  tangent.vals.clear();
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = 0; j < n; ++j) {
      const Real e = T[i * n + j];
      if (e == 0.0) continue;
      tangent.rows.push_back(static_cast<Index>(i));
      tangent.cols.push_back(static_cast<Index>(j));
      tangent.vals.push_back(e);
    }
}

// Advance one HHT step of the NONLINEAR system with a Newton loop on a_{n+1} (task 3.2).
// Predictors carry the old state; each iteration derives u,v from the trial a, assembles
// the material tangent + internal force (+ contact), forms the HHT dynamic residual, and
// solves the effective tangent for the acceleration correction. Returns the iterations
// taken (or -1 on non-convergence within max_iterations).
int step_nonlinear(const Model& model, NonlinDyn& st, const LinearSystem& K,
                   const std::vector<Real>& Md, const std::vector<Real>& Cd,
                   const HhtCoeffs& hc, MotionState& s,
                   const std::vector<Real>& f_new, const std::vector<Real>& f_old,
                   SolverKind kind, std::size_t n, std::size_t n_nodes,
                   const NonlinearControls& ctl) {
  const Real dt = hc.dt, alpha = hc.alpha;
  std::vector<Real> u_p(n), v_p(n);
  for (std::size_t i = 0; i < n; ++i) {
    u_p[i] = s.u[i] + dt * s.v[i] + dt * dt * (0.5 - hc.beta) * s.a[i];
    v_p[i] = s.v[i] + dt * (1.0 - hc.gamma) * s.a[i];
  }
  // Internal force + damping at the START of the step (for the -α old-state HHT term).
  std::vector<Real> fint_old;
  assemble_tangent(model, st, K, s.u, n_nodes, fint_old);
  const std::vector<Real> Cv_old = matvec(Cd, s.v, n);

  MotionState trial = s;  // start from the predictor's a (extrapolate a_{n+1} = a_n)
  newmark_correct(trial, u_p, v_p, s.a, hc, n);
  Real denom = 0.0;
  for (int it = 0; it < ctl.max_iterations; ++it) {
    std::vector<Real> fint_new;
    LinearSystem tangent =
        assemble_tangent(model, st, K, trial.u, n_nodes, fint_new);
    const std::vector<Real> Ma = matvec(Md, trial.a, n);
    const std::vector<Real> Cv_new = matvec(Cd, trial.v, n);
    // HHT residual: (1+α) f_{n+1} - α f_n - M a - (1+α)(f_int_new + C v_new)
    //               + α (f_int_old + C v_old).
    std::vector<Real> r(n);
    for (std::size_t i = 0; i < n; ++i)
      r[i] = (1.0 + alpha) * f_new[i] - alpha * f_old[i] - Ma[i] -
             (1.0 + alpha) * (fint_new[i] + Cv_new[i]) +
             alpha * (fint_old[i] + Cv_old[i]);
    if (it == 0) denom = std::max(norm2(r), ctl.eps);
    if (norm2(r) / denom < ctl.force_tol) {
      s = trial;
      return it + 1;
    }
    fold_effective_tangent(tangent, Md, Cd, hc, n);
    tangent.rhs = r;
    const std::vector<Real> da = solve_reduced(tangent, kind);
    for (std::size_t i = 0; i < n; ++i) trial.a[i] += da[i];
    newmark_correct(trial, u_p, v_p, trial.a, hc, n);
  }
  s = trial;
  return -1;
}

// Damping dense buffer C = alpha M + beta K on the free DOFs.
std::vector<Real> damping_dense(const std::vector<Real>& Md,
                                const std::vector<Real>& Kd, const Damping& d,
                                std::size_t n) {
  std::vector<Real> C(n * n, 0.0);
  if (d.alpha == 0.0 && d.beta == 0.0) return C;
  for (std::size_t k = 0; k < n * n; ++k) C[k] = d.alpha * Md[k] + d.beta * Kd[k];
  return C;
}

// Initial acceleration from rest: M a₀ = f(0) - C v₀ - K u₀ = f(0) (u₀=v₀=0). Solved on
// the mass system (M is PD), so the run starts with a consistent acceleration.
std::vector<Real> initial_acceleration(const LinearSystem& M,
                                       const std::vector<Real>& f0, SolverKind kind) {
  LinearSystem sys = M;
  sys.rhs = f0;
  return solve_reduced(sys, kind);
}

}  // namespace

HhtParams HhtParams::from_alpha(Real alpha) {
  const Real a = std::clamp(alpha, -1.0 / 3.0, 0.0);
  HhtParams p;
  p.alpha = a;
  p.beta = (1.0 - a) * (1.0 - a) / 4.0;
  p.gamma = 0.5 - a;
  return p;
}

std::vector<DirectTimePoint> direct_dynamic(const Model& model, Real dt, Real t_end,
                                            const DirectOptions& opts,
                                            const std::function<Real(Real)>& amplitude,
                                            DirectReport* report) {
  if (dt <= 0.0) throw std::runtime_error("direct_dynamic: dt must be > 0");
  const LinearSystem K = fem::assemble_linear_static(model);
  const LinearSystem M = fem::assemble_mass(model, opts.lumped_mass);
  if (K.n_free != M.n_free)
    throw std::runtime_error("direct_dynamic: K and M free-DOF mismatch");
  const std::size_t n = static_cast<std::size_t>(K.n_free);
  const std::size_t n_nodes =
      K.dof_eq.size() / static_cast<std::size_t>(kDofsPerNode);

  const std::vector<Real> Md = dense_symmetric(M);
  const std::vector<Real> Kd = dense_symmetric(K);
  const std::vector<Real> Cd = damping_dense(Md, Kd, opts.damping, n);
  const std::vector<Real> f_full = fem::external_load_vector(model);
  const std::vector<Real> f_free = reduce_free(K, f_full);

  const HhtCoeffs hc(dt, opts.hht);
  const SolverKind kind =
      opts.forced ? *opts.forced : resolve_solver_kind(model.solver, K.n_free);
  const LinearSystem eff = effective_system(K, Md, Kd, hc, opts.damping, n);

  MotionState s;
  s.u = opts.u0_free.size() == n ? opts.u0_free : std::vector<Real>(n, 0.0);
  s.v = opts.v0_free.size() == n ? opts.v0_free : std::vector<Real>(n, 0.0);
  // M a₀ = f(0) - C v₀ - K u₀ (consistent initial acceleration; reduces to M a₀ = f(0)
  // from rest).
  {
    std::vector<Real> f0 = load_at(f_free, amplitude, 0.0);
    const std::vector<Real> Ku0 = matvec(Kd, s.u, n), Cv0 = matvec(Cd, s.v, n);
    for (std::size_t i = 0; i < n; ++i) f0[i] -= Ku0[i] + Cv0[i];
    s.a = initial_acceleration(M, f0, kind);
  }

  // Nonlinear path (task 3.2): per-step Newton on a_{n+1} with the material/contact
  // tangent folded into the effective dynamic tangent. Built once, advanced per step.
  const bool nl = opts.nonlinear || model.has_nonlinear_material() || model.has_contact();
  NonlinDyn st;
  if (nl) {
    st.mp = fem::make_material_points(model);
    if (model.has_contact()) st.contact = fem::build_contact_pairs(model);
  }

  const std::size_t n_steps = static_cast<std::size_t>(std::llround(t_end / dt));
  std::vector<DirectTimePoint> out;
  out.reserve(n_steps + 1);
  out.push_back(make_frame(K, Md, Kd, s, 0.0, n, n_nodes));

  int total_iters = 0;
  Real t = 0.0;
  for (std::size_t step = 0; step < n_steps; ++step) {
    const std::vector<Real> f_old = load_at(f_free, amplitude, t);
    const std::vector<Real> f_new = load_at(f_free, amplitude, t + dt);
    if (nl) {
      const int iters = step_nonlinear(model, st, K, Md, Cd, hc, s, f_new, f_old,
                                       kind, n, n_nodes, model.controls);
      total_iters += iters > 0 ? iters : model.controls.max_iterations;
      st.mp.commit();  // accept the step's material history
    } else {
      step_linear(eff, Kd, Cd, hc, s, f_new, f_old, kind, n);
    }
    t += dt;
    out.push_back(make_frame(K, Md, Kd, s, t, n, n_nodes));
  }

  if (report) {
    report->steps = static_cast<int>(n_steps);
    report->iterations = total_iters;
    report->nonlinear = nl;
    const Real e0 = out.front().total_energy;
    Real drift = 0.0;
    if (e0 != 0.0)
      for (const DirectTimePoint& tp : out)
        drift = std::max(drift, std::abs(tp.total_energy - e0) / std::abs(e0));
    report->energy_drift = drift;
  }
  return out;
}

}  // namespace cxpp::numerics
