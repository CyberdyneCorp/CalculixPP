#include "calculixpp/numerics/modal_dynamics.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>

#include "numpp/core/creation.hpp"
#include "numpp/core/dtype.hpp"
#include "numpp/core/ndarray.hpp"

// Modal superposition engine (spec: eigensolution / dynamic-analysis). The physics is
// the decoupled SDOF modal equation q̈_k + 2 ζ_k ω_k q̇_k + ω_k² q_k = p_k(t). Transient
// response uses the EXACT piecewise-linear-load recurrence (Nigam-Jennings) so the
// natural period 2π/ω is reproduced with no algorithmic damping; harmonic response uses
// the closed-form complex transfer function 1/(ω_k² - Ω² + 2 i ζ_k ω_k Ω).
namespace cxpp::numerics {
namespace {

// y = M * x on the free DOFs (dense mat-vec), for the base-motion effective load.
std::vector<Real> mass_matvec(const fem::LinearSystem& M, const std::vector<Real>& x) {
  const std::int64_t n = static_cast<std::int64_t>(M.n_free);
  numpp::ndarray A = numpp::zeros({n, n}, numpp::kFloat64);
  for (std::size_t k = 0; k < M.vals.size(); ++k) {
    const std::int64_t r = static_cast<std::int64_t>(M.rows[k]);
    const std::int64_t c = static_cast<std::int64_t>(M.cols[k]);
    A.set_item<double>({r, c}, A.item<double>({r, c}) + M.vals[k]);
  }
  // Symmetrize (the congruence scatter can store one triangle for some entries).
  for (std::int64_t i = 0; i < n; ++i)
    for (std::int64_t j = i + 1; j < n; ++j) {
      const double aij = A.item<double>({i, j}), aji = A.item<double>({j, i});
      const double val = (aij != 0.0 && aji != 0.0) ? aij : aij + aji;
      A.set_item<double>({i, j}, val);
      A.set_item<double>({j, i}, val);
    }
  std::vector<Real> y(static_cast<std::size_t>(n), 0.0);
  for (std::int64_t i = 0; i < n; ++i) {
    double s = 0.0;
    for (std::int64_t j = 0; j < n; ++j)
      s += A.item<double>({i, j}) * x[static_cast<std::size_t>(j)];
    y[static_cast<std::size_t>(i)] = s;
  }
  return y;
}

// Expand a free-DOF vector to a full nodal displacement field through the constraint
// transform (SPC value 0 — a modal response has no prescribed displacement). Mirrors the
// eigensolution shape expansion so the recombined transient field is on the same node
// numbering as a static solve.
std::vector<Vec3> expand_free(const fem::LinearSystem& sys,
                              const std::vector<Real>& q, std::size_t n_nodes) {
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

// One SDOF integrated with the EXACT piecewise-linear-load solution for
//   q̈ + 2ζω q̇ + ω² q = p(τ),   p(τ) = p0 + (p1-p0)/dt · τ   on τ ∈ [0,dt].
// State is (displacement q, velocity qd). The exact solution over the step is the sum
// of a homogeneous part carrying the initial state and a particular part forced by the
// linear load. Both are written in closed form (undamped and under-damped 0≤ζ<1), so the
// natural period 2π/ω is reproduced with zero algorithmic damping — the property the
// analytical SDOF validation relies on. The ω→0 (rigid-body / zero-frequency) mode is the
// free-particle limit q̈ = p, integrated exactly for a linear p.
struct Sdof {
  Real omega{0.0};
  Real zeta{0.0};
  Real q{0.0};
  Real qd{0.0};

  void step(Real dt, Real p0, Real p1) {
    if (omega <= 0.0) {  // rigid-body mode: q̈ = p (double integration of a linear p)
      const Real a0 = p0, a1 = (p1 - p0) / dt;
      const Real q_new = q + qd * dt + a0 * dt * dt / 2.0 + a1 * dt * dt * dt / 6.0;
      const Real qd_new = qd + a0 * dt + a1 * dt * dt / 2.0;
      q = q_new;
      qd = qd_new;
      return;
    }
    const Real w = omega, z = std::min(zeta, 0.999999);
    const Real w2 = w * w;
    const Real wd = w * std::sqrt(std::max(1.0 - z * z, 0.0));  // damped frequency
    const Real e = std::exp(-z * w * dt);
    const Real s = std::sin(wd * dt), c = std::cos(wd * dt);
    const Real slope = (p1 - p0) / dt;

    // Particular solution q_p(τ) = p(τ)/ω² - slope·2ζ/(ω·ω²)  (a line: the quasi-static
    // response of the SDOF to a ramping load). Its value/derivative at the endpoints:
    const Real damp_lag = slope * (2.0 * z) / (w * w2);
    const Real qp0 = p0 / w2 - damp_lag;          // q_p(0)
    const Real qp1 = p1 / w2 - damp_lag;          // q_p(dt)
    const Real qdp = slope / w2;                  // q_p'(τ) (constant)

    // Total state at 0 minus particular gives the transient's initial condition; the
    // transient decays homogeneously and is added back at dt so the particular line
    // shifts the response. h0/hd0 are the homogeneous initial conditions.
    const Real h0 = q - qp0;
    const Real hd0 = qd - qdp;
    // Homogeneous e^{-zω τ}(a cos wd τ + b sin wd τ): a=h0, b=(hd0+zω h0)/wd.
    const Real a = h0;
    const Real b = (hd0 + z * w * h0) / wd;
    const Real h_q = e * (a * c + b * s);
    const Real h_qd =
        e * (-z * w * (a * c + b * s) + wd * (-a * s + b * c));

    q = qp1 + h_q;
    qd = qdp + h_qd;
  }
};

}  // namespace

Real Damping::ratio(std::size_t k, Real omega_k) const {
  if (k < modal_ratios.size()) return modal_ratios[k];
  if (omega_k <= 0.0) return 0.0;
  return 0.5 * (alpha / omega_k + beta * omega_k);
}

std::vector<Real> ModalSystem::project_load(const std::vector<Real>& free_load) const {
  std::vector<Real> p(n_modes, 0.0);
  for (std::size_t k = 0; k < n_modes; ++k) {
    const std::vector<Real>& phi = basis->mode_free[k];
    Real s = 0.0;
    const std::size_t n = std::min(phi.size(), free_load.size());
    for (std::size_t i = 0; i < n; ++i) s += phi[i] * free_load[i];
    p[k] = s;
  }
  return p;
}

ModalSystem project_modal_system(const EigenBasis& basis, const Damping& damping) {
  ModalSystem sys;
  sys.n_modes = basis.modes.size();
  sys.basis = &basis;
  sys.omega.reserve(sys.n_modes);
  sys.lambda.reserve(sys.n_modes);
  sys.zeta.reserve(sys.n_modes);
  for (std::size_t k = 0; k < sys.n_modes; ++k) {
    const Real w = basis.modes[k].omega;
    sys.omega.push_back(w);
    sys.lambda.push_back(w * w);
    sys.zeta.push_back(damping.ratio(k, w));
  }
  return sys;
}

std::vector<ModalTimePoint> modal_dynamic(
    const ModalSystem& sys, const ModalLoad& load, Real dt, Real t_end,
    const std::function<Real(Real)>& amplitude) {
  if (dt <= 0.0) throw std::runtime_error("modal_dynamic: dt must be > 0");
  if (sys.basis == nullptr) throw std::runtime_error("modal_dynamic: null basis");
  const std::vector<Real> p_spatial = sys.project_load(load.pattern);  // φ_kᵀ pattern
  const auto amp = [&](Real t) -> Real { return amplitude ? amplitude(t) : 1.0; };

  std::vector<Sdof> osc(sys.n_modes);
  for (std::size_t k = 0; k < sys.n_modes; ++k) {
    osc[k].omega = sys.omega[k];
    osc[k].zeta = sys.zeta[k];
  }

  const std::size_t n_steps =
      static_cast<std::size_t>(std::llround(t_end / dt));
  const fem::LinearSystem& K = sys.basis->stiffness;
  const std::size_t n_nodes =
      K.dof_eq.size() / static_cast<std::size_t>(kDofsPerNode);

  std::vector<ModalTimePoint> out;
  out.reserve(n_steps + 1);

  auto recombine = [&](Real t) {
    std::vector<Real> q_free(static_cast<std::size_t>(sys.basis->n_free), 0.0);
    for (std::size_t k = 0; k < sys.n_modes; ++k) {
      const std::vector<Real>& phi = sys.basis->mode_free[k];
      const Real qk = osc[k].q;
      for (std::size_t i = 0; i < phi.size(); ++i) q_free[i] += qk * phi[i];
    }
    ModalTimePoint frame;
    frame.time = t;
    frame.displacement = expand_free(K, q_free, n_nodes);
    out.push_back(std::move(frame));
  };

  recombine(0.0);  // rest start
  Real t = 0.0;
  for (std::size_t step = 0; step < n_steps; ++step) {
    const Real a0 = amp(t), a1 = amp(t + dt);
    for (std::size_t k = 0; k < sys.n_modes; ++k)
      osc[k].step(dt, p_spatial[k] * a0, p_spatial[k] * a1);
    t += dt;
    recombine(t);
  }
  return out;
}

HarmonicResponse steady_state_response(const ModalSystem& sys,
                                       const std::vector<Real>& pattern,
                                       Real omega_exc) {
  const std::vector<Real> p = sys.project_load(pattern);  // φ_kᵀ pattern
  HarmonicResponse resp;
  resp.omega = omega_exc;
  resp.frequency = omega_exc / (2.0 * M_PI);
  resp.amplitude.assign(static_cast<std::size_t>(sys.basis->n_free),
                        std::complex<Real>(0.0, 0.0));
  const Real Om = omega_exc;
  for (std::size_t k = 0; k < sys.n_modes; ++k) {
    const Real wk = sys.omega[k], zk = sys.zeta[k];
    // Complex modal amplitude q_k = p_k / (ω_k² - Ω² + 2 i ζ_k ω_k Ω).
    const std::complex<Real> denom(wk * wk - Om * Om, 2.0 * zk * wk * Om);
    if (std::abs(denom) == 0.0) continue;  // undamped exact resonance -> skip (infinite)
    const std::complex<Real> qk = p[k] / denom;
    const std::vector<Real>& phi = sys.basis->mode_free[k];
    for (std::size_t i = 0; i < phi.size(); ++i) resp.amplitude[i] += qk * phi[i];
  }
  return resp;
}

std::vector<HarmonicResponse> steady_state_sweep(const ModalSystem& sys,
                                                 const std::vector<Real>& pattern,
                                                 Real f_lo, Real f_hi,
                                                 std::size_t num) {
  std::vector<HarmonicResponse> out;
  if (num == 0) return out;
  out.reserve(num);
  if (num == 1) {
    out.push_back(steady_state_response(sys, pattern, 2.0 * M_PI * f_lo));
    return out;
  }
  // Logarithmic sweep between f_lo and f_hi (both > 0); linear if f_lo <= 0.
  const bool logs = f_lo > 0.0 && f_hi > 0.0;
  for (std::size_t i = 0; i < num; ++i) {
    const Real frac = static_cast<Real>(i) / static_cast<Real>(num - 1);
    Real f;
    if (logs)
      f = std::exp(std::log(f_lo) + frac * (std::log(f_hi) - std::log(f_lo)));
    else
      f = f_lo + frac * (f_hi - f_lo);
    out.push_back(steady_state_response(sys, pattern, 2.0 * M_PI * f));
  }
  return out;
}

std::vector<Real> base_motion_load(const EigenBasis& basis,
                                   const fem::LinearSystem& M, int dir) {
  const std::int64_t n = static_cast<std::int64_t>(basis.n_free);
  const std::size_t n_nodes =
      basis.stiffness.dof_eq.size() / static_cast<std::size_t>(kDofsPerNode);
  // Rigid-body influence vector r for direction `dir` on the free DOFs (unit
  // translation at every node, projected through the constraint transform).
  std::vector<Real> r(static_cast<std::size_t>(n), 0.0);
  const fem::ConstraintTransform& tf = basis.stiffness.transform;
  for (std::size_t ni = 0; ni < n_nodes; ++ni) {
    const std::size_t g = ni * kDofsPerNode + static_cast<std::size_t>(dir);
    if (g >= tf.expansion.size()) continue;
    for (const fem::DofTerm& t : tf.expansion[g].terms)
      r[static_cast<std::size_t>(t.eq)] += t.coeff;
  }
  // Effective load f_eff = -M r (so p_k = φ_kᵀ f_eff = -Γ_k).
  std::vector<Real> Mr = mass_matvec(M, r);
  for (Real& v : Mr) v = -v;
  return Mr;
}

}  // namespace cxpp::numerics
