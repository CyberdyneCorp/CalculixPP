#include "calculixpp/numerics/heat_transfer.hpp"

#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "calculixpp/fem/cavity_radiation.hpp"
#include "calculixpp/fem/element.hpp"
#include "calculixpp/fem/loads.hpp"
#include "calculixpp/fem/thermal.hpp"
#include "calculixpp/numerics/linear_static.hpp"
#include "calculixpp/numerics/nonlinear_static.hpp"

namespace cxpp::numerics {
namespace {

using fem::FaceSurface;
using fem::LinearSystem;

// A full-node sparse operator accumulated as a COO map keyed by (row,col) over the
// FULL node numbering, plus a full-node rhs. Reduced to the free system (Dirichlet
// elimination) just before each solve so radiation/transient drivers can rebuild the
// operator cheaply per iteration without re-numbering DOFs.
struct FullOperator {
  std::unordered_map<std::int64_t, Real> a;  // (row*N + col) -> value
  std::vector<Real> b;                        // rhs per node, size N
  Index N{};

  void add(Index r, Index c, Real v) {
    if (v != 0.0) a[static_cast<std::int64_t>(r) * N + c] += v;
  }
};

FullOperator make_full_operator(Index n_nodes) {
  FullOperator op;
  op.N = n_nodes;
  op.b.assign(static_cast<std::size_t>(n_nodes), 0.0);
  return op;
}

// Reduce a full-node operator A x = b into the free/free LinearSystem using the
// scalar thermal DOF map already built in `sys` (dof_eq / prescribed / n_free),
// moving prescribed columns to the rhs. Rebuilds sys.rows/cols/vals/rhs.
void reduce_operator(const FullOperator& op, LinearSystem& sys) {
  const Index nf = sys.n_free;
  std::unordered_map<std::int64_t, Real> kmap;
  sys.rhs.assign(static_cast<std::size_t>(nf), 0.0);
  for (const auto& [key, val] : op.a) {
    const Index r = static_cast<Index>(key / op.N);
    const Index c = static_cast<Index>(key % op.N);
    const Index eqr = sys.dof_eq[static_cast<std::size_t>(r)];
    if (eqr < 0) continue;  // prescribed row eliminated
    const Index eqc = sys.dof_eq[static_cast<std::size_t>(c)];
    if (eqc >= 0)
      kmap[static_cast<std::int64_t>(eqr) * nf + eqc] += val;
    else
      sys.rhs[static_cast<std::size_t>(eqr)] -= val * sys.prescribed[static_cast<std::size_t>(c)];
  }
  for (Index ni = 0; ni < op.N; ++ni) {
    const Index eq = sys.dof_eq[static_cast<std::size_t>(ni)];
    if (eq >= 0) sys.rhs[static_cast<std::size_t>(eq)] += op.b[static_cast<std::size_t>(ni)];
  }
  sys.rows.clear();
  sys.cols.clear();
  sys.vals.clear();
  sys.rows.reserve(kmap.size());
  sys.cols.reserve(kmap.size());
  sys.vals.reserve(kmap.size());
  for (const auto& [k, v] : kmap) {
    sys.rows.push_back(static_cast<Index>(k / nf));
    sys.cols.push_back(static_cast<Index>(k % nf));
    sys.vals.push_back(v);
  }
}

// Scatter the free-DOF solution back onto the full nodal temperature field.
std::vector<Real> expand_solution(const LinearSystem& sys,
                                  const std::vector<Real>& tf) {
  std::vector<Real> temp(sys.dof_eq.size(), 0.0);
  for (std::size_t ni = 0; ni < sys.dof_eq.size(); ++ni) {
    const Index eq = sys.dof_eq[ni];
    temp[ni] = eq >= 0 ? tf[static_cast<std::size_t>(eq)] : sys.prescribed[ni];
  }
  return temp;
}

// Add the linearized surface-to-ambient radiation of every *RADIATE face to the
// operator around the current temperature `temp` (full residual form): for each
// face node a, with T_abs = T - absolute_zero and c = eps*sigma:
//   flux out  q_a = ∫ N_a c (T_abs^4 - Tamb^4) dA
//   tangent   dq_a/dT_b = ∫ N_a c 4 T_abs^3 N_b dA   (via the face mass matrix)
// We move the tangent to the LHS and (tangent*T - residual) to the RHS so the
// solve yields the updated temperature directly (Newton on A(T) T = b): the net
// contribution to b at node a is  Σ_b tangent_ab T_b - q_a. Radiation is inert when
// sigma == 0. Returns the max |residual| for a convergence check.
Real add_radiation(const Model& model, const std::vector<Real>& temp, Real lambda,
                   FullOperator& op) {
  const Real sigma = model.physical.sigma;
  if (sigma == 0.0 || model.radiates.empty()) return 0.0;
  const Real az = model.physical.absolute_zero;
  const Mesh& mesh = model.mesh;
  Real max_res = 0.0;

  for (const Radiate& rd : model.radiates) {
    if (rd.cavity) continue;  // gray-body cavity exchange is added separately
    const Index ei = mesh.element_index(rd.elem_id);
    if (ei < 0) throw std::runtime_error("*RADIATE references unknown element");
    const FaceSurface fs = fem::face_surface_integrals(mesh, ei, rd.face);
    const int nf = static_cast<int>(fs.gnode.size());
    const Real c = rd.emissivity * sigma;
    const Real tamb = model.amplitude_factor(rd.amplitude, lambda) * rd.ambient_temp - az;
    const Real tamb4 = tamb * tamb * tamb * tamb;

    for (int a = 0; a < nf; ++a) {
      const std::size_t na = static_cast<std::size_t>(fs.gnode[static_cast<std::size_t>(a)]);
      Real qa = 0.0, tan_dot_t = 0.0;
      for (int b = 0; b < nf; ++b) {
        const std::size_t nb = static_cast<std::size_t>(fs.gnode[static_cast<std::size_t>(b)]);
        const Real m = fs.mass[static_cast<std::size_t>(a) * static_cast<std::size_t>(nf) +
                               static_cast<std::size_t>(b)];
        const Real tb = temp[nb] - az;              // T_abs at node b
        // q_a uses the field value; approximate T_abs^4 by node-b interpolation of
        // T_abs then raise to 4th at the node (nodal collocation of the quartic).
        const Real tb4 = tb * tb * tb * tb;
        qa += c * m * (tb4 - tamb4);
        const Real dq = c * 4.0 * tb * tb * tb * m;  // dq_a/dT_b
        tan_dot_t += dq * temp[nb];
        op.add(fs.gnode[static_cast<std::size_t>(a)], fs.gnode[static_cast<std::size_t>(b)], dq);
      }
      op.b[na] += tan_dot_t - qa;  // Newton move: (K_t T - R) to the rhs
      if (std::abs(qa) > max_res) max_res = std::abs(qa);
    }
  }
  return max_res;
}

// Add the gray-body cavity (surface-to-surface) radiation exchange to the operator,
// Newton-linearized around the current temperature `temp`. Each cavity patch i has an
// area-consistent mean absolute temperature T_i = (Σ_a load_a T_a)/A_i; the radiosity
// solve gives the net heat LOSS Q_i(T) leaving the patch and its tangent dQ_i/dT_k
// (see fem::cavity_heat_flow). The loss is distributed to the patch nodes through the
// face load shape (q per area = Q_i/A_i, node a gets (Q_i/A_i) load_a), and the tangent
// couples patch temperatures via dT_k/dT_b = load_b/A_k. The residual added to node a
// of patch i is  R_a = (Q_i/A_i) load_a  (a heat sink), moved to the rhs as
// (Σ tangent*T - R). Returns the max |Q_i| for the convergence check. Cavity geometry
// (view factors) is fixed, so `cav` is built once by the caller.
Real add_cavity_radiation(const Model& model, const fem::Cavity& cav,
                          const std::vector<Real>& temp, FullOperator& op) {
  const Real sigma = model.physical.sigma;
  if (sigma == 0.0 || cav.n == 0) return 0.0;
  const Real az = model.physical.absolute_zero;
  const int n = cav.n;

  // Patch mean absolute temperatures T_i = (Σ_a load_a (temp_a - az)) / A_i.
  std::vector<Real> tabs(static_cast<std::size_t>(n), 0.0);
  for (int i = 0; i < n; ++i) {
    const fem::CavityPatch& p = cav.patches[static_cast<std::size_t>(i)];
    Real acc = 0.0;
    for (std::size_t a = 0; a < p.surf.gnode.size(); ++a)
      acc += p.surf.load[a] * (temp[static_cast<std::size_t>(p.surf.gnode[a])] - az);
    tabs[static_cast<std::size_t>(i)] = p.area > 0.0 ? acc / p.area : 0.0;
  }

  std::vector<Real> Q, dQdT;
  fem::cavity_heat_flow(cav, tabs, sigma, Q, dQdT);

  Real max_res = 0.0;
  for (int i = 0; i < n; ++i) {
    const fem::CavityPatch& p = cav.patches[static_cast<std::size_t>(i)];
    if (p.area <= 0.0) continue;
    const Real qi = Q[static_cast<std::size_t>(i)] / p.area;  // net loss per unit area
    if (std::abs(Q[static_cast<std::size_t>(i)]) > max_res) max_res = std::abs(Q[static_cast<std::size_t>(i)]);

    for (std::size_t a = 0; a < p.surf.gnode.size(); ++a) {
      const Index na = p.surf.gnode[a];
      const Real Ra = qi * p.surf.load[a];  // residual (heat sink) at node a
      Real tan_dot_t = 0.0;
      // Tangent: dRa/dT_b(patch k) = (load_a/A_i) dQ_i/dT_k (load_b/A_k). Assemble the
      // full coupling and simultaneously accumulate (tangent . T) for the Newton rhs.
      for (int k = 0; k < n; ++k) {
        const fem::CavityPatch& pk = cav.patches[static_cast<std::size_t>(k)];
        if (pk.area <= 0.0) continue;
        const Real dq = dQdT[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) +
                             static_cast<std::size_t>(k)];
        const Real coef = (p.surf.load[a] / p.area) * dq / pk.area;
        for (std::size_t b = 0; b < pk.surf.gnode.size(); ++b) {
          const Index nb = pk.surf.gnode[b];
          const Real kab = coef * pk.surf.load[b];
          if (kab == 0.0) continue;
          op.add(na, nb, kab);
          tan_dot_t += kab * temp[static_cast<std::size_t>(nb)];
        }
      }
      op.b[static_cast<std::size_t>(na)] += tan_dot_t - Ra;  // (K_t T - R) to the rhs
    }
  }
  return max_res;
}

// Seed the full conduction+film operator from the shared FEM assembly at step
// fraction `lambda` into `op`.
void seed_conduction(const fem::FullThermalSystem& fts, FullOperator& op) {
  for (std::size_t i = 0; i < fts.k_vals.size(); ++i)
    op.add(fts.k_rows[i], fts.k_cols[i], fts.k_vals[i]);
  for (Index ni = 0; ni < op.N; ++ni)
    op.b[static_cast<std::size_t>(ni)] += fts.flux[static_cast<std::size_t>(ni)];
}

// Heat-flux reaction RFL = Kt_conduction * T at every node (CalculiX "heat
// generation"). This is the pure-conduction internal flux: film, radiation, and the
// capacitance term are NOT included, matching the CalculiX RFL convention validated
// against oneel20cf (conduction) and oneel20fi (film — RFL is still conduction-only).
std::vector<Real> flux_reactions(const Model& model,
                                 const std::vector<Real>& temperature) {
  const Mesh& mesh = model.mesh;
  const std::vector<Real> cond = model.has_temp_dependent_thermal()
                                     ? model.element_conductivity_at(temperature)
                                     : model.element_conductivity();
  std::vector<Real> rfl(mesh.num_nodes(), 0.0);

  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    const Element& elem = mesh.elements()[e];
    const int n = nodes_per_element(elem.type);
    std::vector<Vec3> coords(static_cast<std::size_t>(n));
    std::vector<Index> nidx(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
      const Index ni = mesh.node_index(elem.nodes[static_cast<std::size_t>(i)]);
      nidx[static_cast<std::size_t>(i)] = ni;
      coords[static_cast<std::size_t>(i)] = mesh.nodes()[static_cast<std::size_t>(ni)].x;
    }
    const std::vector<Real> Kt = fem::element_conduction(elem.type, coords, cond[e]);
    for (int a = 0; a < n; ++a) {
      Real v = 0.0;
      for (int b = 0; b < n; ++b)
        v += Kt[static_cast<std::size_t>(a) * static_cast<std::size_t>(n) +
                static_cast<std::size_t>(b)] *
             temperature[static_cast<std::size_t>(nidx[static_cast<std::size_t>(b)])];
      rfl[static_cast<std::size_t>(nidx[static_cast<std::size_t>(a)])] += v;
    }
  }
  return rfl;
}

// Iterate the nonlinear (radiation + temperature-dependent conductivity) solve to
// convergence. `build_base` rebuilds the linear operator (conduction + film [+ C/dt])
// from the CURRENT temperature field, so a temperature-dependent k(T)/c(T) is
// re-evaluated (Picard) each iteration; radiation is re-linearized on top of it. With
// constant properties and no radiation the base is fixed and this converges in one
// iteration (the base is rebuilt but identical), so the constant path is unchanged.
// `temp` is the current field, updated in place. Shared by the steady solve and each
// transient backward-Euler step.
template <typename BuildBase>
void nonlinear_solve(const Model& model, const fem::Cavity& cav, LinearSystem& sys,
                     BuildBase&& build_base, SolverKind kind, Real lambda,
                     std::vector<Real>& temp) {
  const bool has_rad = model.physical.sigma != 0.0 &&
                       (!model.radiates.empty() || cav.n > 0);
  const bool iterate = has_rad || model.has_temp_dependent_thermal();
  const int max_it = iterate ? 40 : 1;
  for (int it = 0; it < max_it; ++it) {
    FullOperator op = build_base(temp);
    const Real res = add_radiation(model, temp, lambda, op) +
                     add_cavity_radiation(model, cav, temp, op);
    reduce_operator(op, sys);
    const std::vector<Real> next = expand_solution(sys, solve_reduced(sys, kind));
    Real dmax = 0.0;
    for (std::size_t i = 0; i < next.size(); ++i)
      dmax = std::max(dmax, std::abs(next[i] - temp[i]));
    temp = next;
    if (!iterate || (res < 1e-9 && dmax < 1e-9)) break;
  }
}

ThermalFields finalize(const Model& model, std::vector<Real> temperature) {
  ThermalFields res;
  res.flux_reaction = flux_reactions(model, temperature);
  fem::recover_heat_flux(model, temperature, res);  // HFL (GP + nodal)
  res.temperature = std::move(temperature);
  return res;
}

// The initial temperature field for a transient step: uniform default + per-node
// *INITIAL CONDITIONS overrides, with prescribed-temperature nodes pinned to their
// BC value so the first backward-Euler step starts consistent.
std::vector<Real> initial_field(const Model& model, const LinearSystem& sys) {
  const Mesh& mesh = model.mesh;
  std::vector<Real> temp(mesh.num_nodes(), model.initial_temperature);
  for (const auto& [node_id, val] : model.initial_temperature_by_node) {
    const Index ni = mesh.node_index(node_id);
    if (ni >= 0) temp[static_cast<std::size_t>(ni)] = val;
  }
  for (std::size_t ni = 0; ni < sys.dof_eq.size(); ++ni)
    if (sys.dof_eq[ni] < 0) temp[ni] = sys.prescribed[ni];
  return temp;
}

// Build one backward-Euler step operator: conduction + film + C/h on the LHS and
// (C/h) T_n added to the rhs, at step fraction `lambda` and step size `h`. Temperature-
// dependent conductivity/capacitance are evaluated at `temp_cond` (the current Picard
// iterate); the capacitance rhs (C/h) T_n uses the previous-step field `temp_prev`.
FullOperator step_operator(const Model& model, Index N, Real lambda, Real h,
                           const std::vector<Real>& temp_prev,
                           const std::vector<Real>& temp_cond) {
  const fem::FullThermalSystem fts = fem::assemble_full_thermal(model, lambda, temp_cond);
  FullOperator base = make_full_operator(N);
  seed_conduction(fts, base);
  for (std::size_t i = 0; i < fts.c_vals.size(); ++i) {
    const Real cij = fts.c_vals[i] / h;
    base.add(fts.c_rows[i], fts.c_cols[i], cij);
    base.b[static_cast<std::size_t>(fts.c_rows[i])] +=
        cij * temp_prev[static_cast<std::size_t>(fts.c_cols[i])];
  }
  return base;
}

// Transient backward-Euler: (C/dt + Kt_film) T_{n+1} = (C/dt) T_n + q(+radiation),
// stepped from the initial field over the step period (increment.total) with the
// step's initial increment. Radiation is re-linearized inside each step.
ThermalFields solve_transient(const Model& model, const fem::Cavity& cav,
                              LinearSystem& sys, SolverKind kind) {
  const Index N = static_cast<Index>(model.mesh.num_nodes());
  const Real period = model.increment.total > 0.0 ? model.increment.total : 1.0;
  Real dt = model.increment.initial > 0.0 ? model.increment.initial : period;
  if (dt > period) dt = period;

  std::vector<Real> temp = initial_field(model, sys);
  Real t = 0.0;
  const int max_steps = 1000000;
  for (int step = 0; step < max_steps && t < period - 1e-12; ++step) {
    const Real h = std::min(dt, period - t);
    const Real lambda = (t + h) / period;
    const std::vector<Real> temp_prev = temp;  // T_n for the (C/h) T_n capacitance rhs
    nonlinear_solve(
        model, cav, sys,
        [&](const std::vector<Real>& cur) {
          return step_operator(model, N, lambda, h, temp_prev, cur);
        },
        kind, lambda, temp);
    t += h;
  }
  return finalize(model, temp);
}

}  // namespace

ThermalFields solve_heat_transfer(const Model& model,
                                  std::optional<SolverKind> forced) {
  LinearSystem sys;
  fem::build_thermal_dof_map(model, sys);
  const SolverKind kind =
      forced ? *forced : resolve_solver_kind(model.solver, sys.n_free);

  // Cavity geometry (view factors) is fixed for the step, so build it once and share
  // it across every Newton/transient iteration. Empty (cav.n == 0) when no *RADIATE
  // ...,CR face is present, in which case the cavity path is inert.
  const fem::Cavity cav = fem::build_cavity(model);

  if (model.procedure == Procedure::HeatTransferTransient)
    return solve_transient(model, cav, sys, kind);

  // Steady state: conduction + film at full magnitude (lambda = 1), radiation and
  // temperature-dependent k(T) iterated. The base is rebuilt from the current field
  // each iteration (constant properties -> identical base -> one iteration).
  const Index N = static_cast<Index>(model.mesh.num_nodes());
  std::vector<Real> temp = expand_solution(
      sys, std::vector<Real>(static_cast<std::size_t>(sys.n_free), 0.0));
  nonlinear_solve(
      model, cav, sys,
      [&](const std::vector<Real>& cur) {
        const fem::FullThermalSystem fts = fem::assemble_full_thermal(model, 1.0, cur);
        FullOperator base = make_full_operator(N);
        seed_conduction(fts, base);
        return base;
      },
      kind, 1.0, temp);
  return finalize(model, std::move(temp));
}

namespace {

// Build the mechanical sub-model for the displacement solve of a coupled step: the
// same model with the thermal procedure swapped to Static and its applied_temperature
// field set to the current nodal temperatures `temp` (so the thermal strain
// eps_th = alpha (T - Tref) drives the mechanical solve). `plastic_heat` (per element,
// optional) is the plastic-dissipation heat source fed back into the thermal solve; it
// is not used here.
Model mechanical_submodel(const Model& model, const std::vector<Real>& temp) {
  Model mech = model;
  mech.procedure = Procedure::Static;
  mech.applied_temperature.clear();
  const Mesh& mesh = model.mesh;
  for (std::size_t ni = 0; ni < mesh.num_nodes(); ++ni)
    mech.applied_temperature[mesh.nodes()[ni].id] = temp[ni];
  return mech;
}

// Run the displacement solve of a coupled step at temperature `temp`. Uses the
// nonlinear driver when the model carries a nonlinear material (so plasticity — the
// source of the two-way heat — is integrated), otherwise the linear solve. Fills
// `eqplastic` (per element) when non-null so the caller can compute the dissipation
// heat. The two paths agree exactly for a linear-elastic thermal-stress model (the
// nonlinear driver reproduces the linear solve), so the scheme's decoupled result is
// unchanged.
StaticFields mechanical_solve(const Model& mech, std::optional<SolverKind> forced,
                              std::vector<Real>* eqplastic) {
  if (mech.has_nonlinear_material()) {
    NonlinearOptions opts;
    opts.forced = forced;
    return solve_nonlinear_static(mech, opts, nullptr, eqplastic);
  }
  if (eqplastic) eqplastic->assign(mech.mesh.num_elements(), 0.0);
  return solve_linear_static(mech, forced);
}

// Inject the per-element plastic-dissipation heat `q_elem` (watts per element) into a
// thermal sub-model as concentrated nodal fluxes (*CFLUX-equivalent): each element's
// heat is split equally over its nodes (a lumped body source, energy-conserving —
// Sum over nodes == Q_e). Returns a model copy whose cfluxes carry the extra source;
// with q_elem all-zero it is the input model unchanged (one-way coupling).
Model with_dissipation_source(const Model& model, const std::vector<Real>& q_elem) {
  bool any = false;
  for (const Real q : q_elem)
    if (q != 0.0) { any = true; break; }
  if (!any) return model;
  Model th = model;
  const Mesh& mesh = model.mesh;
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    if (q_elem[e] == 0.0) continue;
    const Element& el = mesh.elements()[e];
    const int nn = nodes_per_element(el.type);
    const Real per_node = q_elem[e] / static_cast<Real>(nn);
    for (int i = 0; i < nn; ++i)
      th.cfluxes.push_back(Cflux{el.nodes[static_cast<std::size_t>(i)], per_node, ""});
  }
  return th;
}

// Max absolute difference between two nodal-temperature fields (coupled convergence).
Real field_change(const std::vector<Real>& a, const std::vector<Real>& b) {
  Real m = 0.0;
  for (std::size_t i = 0; i < a.size() && i < b.size(); ++i)
    m = std::max(m, std::abs(a[i] - b[i]));
  return m;
}

// Gauss-Seidel STAGGERED coupled solve (task 4.2). Alternate a thermal solve and a
// mechanical solve, feeding the mechanical plastic-dissipation heat back into the next
// thermal solve, until the temperature field stops changing. Converges in ONE pass
// when the coupling is one-way (taylor_quinney == 0 or no plasticity): the first
// thermal solve has no dissipation source, and the second iteration's thermal field is
// identical, so the loop exits after one mechanical solve — byte-for-byte the
// historical sequential result. `report` (optional) records the outer iteration count.
CoupledFields solve_coupled_staggered(const Model& model,
                                      std::optional<SolverKind> forced,
                                      int* outer_iters) {
  CoupledFields out;
  std::vector<Real> q_elem(model.mesh.num_elements(), 0.0);
  std::vector<Real> prev_temp;
  const int max_outer = 30;
  int it = 0;
  for (; it < max_outer; ++it) {
    const Model th = with_dissipation_source(model, q_elem);
    out.thermal = solve_heat_transfer(th, forced);

    const Model mech = mechanical_submodel(model, out.thermal.temperature);
    std::vector<Real> eqplastic;
    out.mechanical = mechanical_solve(mech, forced, &eqplastic);
    q_elem = model.plastic_dissipation_heat(eqplastic);

    const bool converged = !prev_temp.empty() &&
                           field_change(out.thermal.temperature, prev_temp) < 1e-9;
    prev_temp = out.thermal.temperature;
    // One-way: q_elem stays all-zero, so the next thermal solve is identical -> the
    // second iteration converges immediately. Two-way: iterate until the fed-back heat
    // no longer moves the temperature field.
    bool has_source = false;
    for (const Real q : q_elem)
      if (q != 0.0) { has_source = true; break; }
    if ((!has_source && it >= 0) || converged) { ++it; break; }
  }
  if (outer_iters) *outer_iters = it;
  return out;
}

// Sparse COO accumulator for the combined 4-DOF/node coupled system, keyed by the
// combined free-DOF numbering (mechanical [0,nf_m) then thermal [nf_m, nf_m+nf_t)).
struct CoupledCoo {
  Index N{};
  std::unordered_map<std::int64_t, Real> a;
  std::vector<Real> rhs;
  void add(Index r, Index c, Real v) {
    if (v != 0.0) a[static_cast<std::int64_t>(r) * N + c] += v;
  }
};

// Assemble the K_uT thermal-strain coupling block into the combined system: for each
// element with *EXPANSION, -C_e maps nodal (T_j - Tref) to mechanical DOF a. The row is
// reduced through the mechanical constraint transform, the column through the thermal
// Dirichlet map; a prescribed temperature and the -Tref constant land on the rhs.
void add_thermal_coupling(const Model& mech, const fem::LinearSystem& msys,
                          const fem::LinearSystem& tsys, Index nf_m, CoupledCoo& co) {
  const Mesh& mesh = mech.mesh;
  const std::vector<ElasticIso> elastic = mech.element_elastic();
  const std::vector<std::optional<Expansion>> expansion = mech.element_expansion();
  const std::vector<bool> active = mech.element_active_mask();
  for (std::size_t e = 0; e < mesh.num_elements(); ++e) {
    if (!active[e] || !expansion[e] || expansion[e]->empty()) continue;
    const Element& el = mesh.elements()[e];
    const int n = nodes_per_element(el.type);
    std::vector<Vec3> coords(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i)
      coords[static_cast<std::size_t>(i)] =
          mesh.nodes()[static_cast<std::size_t>(
              mesh.node_index(el.nodes[static_cast<std::size_t>(i)]))].x;
    // K_Tu == 0 here, so the coupling block does not feed back into T; the constant
    // (first) alpha is exact for a constant-alpha deck and the mechanical re-solve
    // below recovers the exact alpha(T) stress otherwise.
    const Real alpha =
        expansion[e]->alpha.value.empty() ? 0.0 : expansion[e]->alpha.value[0];
    const fem::D6 D = fem::elastic_iso_D(elastic[e]);
    const std::vector<Real> C = fem::element_thermal_coupling(el.type, coords, D, alpha);
    const Real tref = expansion[e]->t_ref;
    for (int a = 0; a < n * kDofsPerNode; ++a) {
      const std::size_t g_a =
          static_cast<std::size_t>(mesh.node_index(
              el.nodes[static_cast<std::size_t>(a / kDofsPerNode)])) *
              kDofsPerNode + static_cast<std::size_t>(a % kDofsPerNode);
      const fem::DofExpansion& ex = msys.transform.expansion[g_a];
      if (ex.terms.empty()) continue;  // prescribed/eliminated mechanical row
      for (int jn = 0; jn < n; ++jn) {
        const Real cval = C[static_cast<std::size_t>(a) * static_cast<std::size_t>(n) +
                            static_cast<std::size_t>(jn)];
        if (cval == 0.0) continue;
        const Index nj = mesh.node_index(el.nodes[static_cast<std::size_t>(jn)]);
        const Index eqt = tsys.dof_eq[static_cast<std::size_t>(nj)];
        for (const fem::DofTerm& rt : ex.terms) {
          if (eqt >= 0)
            co.add(rt.eq, nf_m + eqt, -rt.coeff * cval);
          else
            co.rhs[static_cast<std::size_t>(rt.eq)] +=
                rt.coeff * cval * tsys.prescribed[static_cast<std::size_t>(nj)];
          co.rhs[static_cast<std::size_t>(rt.eq)] += -rt.coeff * cval * tref;
        }
      }
    }
  }
}

// Assemble and solve the MONOLITHIC 4-DOF/node (u,v,w,T) coupled system in ONE linear
// solve (task 4.1), for the given per-element plastic-dissipation heat `q_elem` (the
// mechanical->thermal source; all-zero for the one-way elastic case). The combined
// reduced system stacks the mechanical free block [0,nf_m) and the thermal free block
// [nf_m, nf_m+nf_t):
//     [ K_uu  K_uT ] [u]   [ f_ext ]
//     [  0    K_TT ] [T] = [  q     ]
// K_uu / K_TT are the already-reduced mechanical / thermal operators; K_uT is the
// thermal-strain coupling C_e (moved to the LHS with a minus sign), reduced through the
// mechanical constraint transform (rows) and the thermal Dirichlet map (cols) — a
// prescribed temperature or the -Tref constant lands on the mechanical rhs. K_Tu == 0
// here (dissipation enters only through `q` on the rhs), so the system is block-
// triangular and its solution equals the sequential thermal-then-mechanical solve
// EXACTLY — the required validation — while being one assembled 4-DOF operator.
CoupledFields solve_monolithic_once(const Model& model,
                                    std::optional<SolverKind> forced,
                                    const std::vector<Real>& q_elem,
                                    std::vector<Real>* eqplastic) {
  const Mesh& mesh = model.mesh;
  // Thermal reduced system (K_TT, thermal Dirichlet map), with the dissipation source.
  const Model th = with_dissipation_source(model, q_elem);
  fem::LinearSystem tsys;
  fem::build_thermal_dof_map(th, tsys);
  const fem::LinearSystem tred = fem::assemble_conduction(th);
  // Mechanical reduced system (K_uu, constraint transform) at zero temperature — the
  // thermal-strain load is supplied instead through the K_uT coupling block below.
  Model mech = model;
  mech.procedure = Procedure::Static;
  mech.applied_temperature.clear();
  const fem::LinearSystem msys = fem::assemble_linear_static(mech);

  const Index nf_m = msys.n_free, nf_t = tsys.n_free;
  const Index N = nf_m + nf_t;
  CoupledCoo co;
  co.N = N;
  co.rhs.assign(static_cast<std::size_t>(N), 0.0);
  // Diagonal blocks + their reduced rhs (external load / thermal flux).
  for (std::size_t i = 0; i < msys.vals.size(); ++i)
    co.add(msys.rows[i], msys.cols[i], msys.vals[i]);
  for (Index i = 0; i < nf_m; ++i) co.rhs[static_cast<std::size_t>(i)] = msys.rhs[static_cast<std::size_t>(i)];
  for (std::size_t i = 0; i < tred.vals.size(); ++i)
    co.add(nf_m + tred.rows[i], nf_m + tred.cols[i], tred.vals[i]);
  for (Index i = 0; i < nf_t; ++i)
    co.rhs[static_cast<std::size_t>(nf_m + i)] = tred.rhs[static_cast<std::size_t>(i)];

  add_thermal_coupling(mech, msys, tsys, nf_m, co);

  // Solve the combined system.
  fem::LinearSystem comb;
  comb.n_free = N;
  comb.rhs = co.rhs;
  comb.rows.reserve(co.a.size());
  comb.cols.reserve(co.a.size());
  comb.vals.reserve(co.a.size());
  for (const auto& [k, v] : co.a) {
    comb.rows.push_back(static_cast<Index>(k / N));
    comb.cols.push_back(static_cast<Index>(k % N));
    comb.vals.push_back(v);
  }
  const SolverKind kind = forced ? *forced : resolve_solver_kind(model.solver, N);
  const std::vector<Real> x = solve_reduced(comb, kind);

  // Split the solution: thermal free DOFs -> temperature field; then expand + recover.
  std::vector<Real> tf(static_cast<std::size_t>(nf_t), 0.0);
  for (Index i = 0; i < nf_t; ++i) tf[static_cast<std::size_t>(i)] = x[static_cast<std::size_t>(nf_m + i)];
  CoupledFields out;
  out.thermal = finalize(th, expand_solution(tsys, tf));

  // Rebuild the mechanical fields from the solved temperature so stress recovery uses
  // sigma = D (eps_mech - eps_th): run the standard mechanical solve at that T. Because
  // the K_uu block and the coupling reproduce K_uu u = f_ext + f_th(T), the displacement
  // from that solve equals the monolithic u to solver tolerance (validated).
  const Model mech_at_T = mechanical_submodel(model, out.thermal.temperature);
  out.mechanical = mechanical_solve(mech_at_T, forced, eqplastic);
  return out;
}

// MONOLITHIC coupled solve (task 4.1): assemble+solve the 4-DOF/node system once (no
// two-way term), or wrap the single solve in an outer fixed point when plastic-
// dissipation heating (the mechanical->thermal K_Tu source) is present, so the joint
// (T,u) state converges. The thermal field carries no MPCs, so the merged SPC/coupling
// reduction always applies.
CoupledFields solve_coupled_monolithic(const Model& model,
                                       std::optional<SolverKind> forced,
                                       int* outer_iters) {
  CoupledFields out;
  std::vector<Real> q_elem(model.mesh.num_elements(), 0.0);
  std::vector<Real> prev_temp;
  const int max_outer = 30;
  int it = 0;
  for (; it < max_outer; ++it) {
    std::vector<Real> eqplastic;
    out = solve_monolithic_once(model, forced, q_elem, &eqplastic);
    q_elem = model.plastic_dissipation_heat(eqplastic);
    bool has_source = false;
    for (const Real q : q_elem)
      if (q != 0.0) { has_source = true; break; }
    const bool converged = !prev_temp.empty() &&
                           field_change(out.thermal.temperature, prev_temp) < 1e-9;
    prev_temp = out.thermal.temperature;
    if (!has_source || converged) { ++it; break; }
  }
  if (outer_iters) *outer_iters = it;
  return out;
}

}  // namespace

CoupledFields solve_coupled(const Model& model, std::optional<SolverKind> forced) {
  int outer = 0;
  return model.coupled_scheme == CoupledScheme::Monolithic
             ? solve_coupled_monolithic(model, forced, &outer)
             : solve_coupled_staggered(model, forced, &outer);
}

CoupledReport solve_coupled_reported(const Model& model,
                                     std::optional<SolverKind> forced) {
  CoupledReport rep;
  rep.fields = model.coupled_scheme == CoupledScheme::Monolithic
                   ? solve_coupled_monolithic(model, forced, &rep.outer_iterations)
                   : solve_coupled_staggered(model, forced, &rep.outer_iterations);
  return rep;
}

}  // namespace cxpp::numerics
