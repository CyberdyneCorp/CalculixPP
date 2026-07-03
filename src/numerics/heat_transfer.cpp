#include "calculixpp/numerics/heat_transfer.hpp"

#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "calculixpp/fem/element.hpp"
#include "calculixpp/fem/loads.hpp"
#include "calculixpp/fem/thermal.hpp"
#include "calculixpp/numerics/linear_static.hpp"

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
  const std::vector<Real> cond = model.element_conductivity();
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

// Iterate the nonlinear (radiation) steady solve to convergence around a fixed
// conduction+film operator. `base` is seeded once; radiation is re-linearized each
// iteration. Converges in one iteration when there is no radiation. `base` is the
// linear operator (conduction + film [+ C/dt for a transient step]); `temp` is the
// current field, updated in place. Shared by the steady solve and each transient
// backward-Euler step.
void nonlinear_solve(const Model& model, LinearSystem& sys, const FullOperator& base,
                     SolverKind kind, Real lambda, std::vector<Real>& temp) {
  const bool has_rad = model.physical.sigma != 0.0 && !model.radiates.empty();
  const int max_it = has_rad ? 40 : 1;
  for (int it = 0; it < max_it; ++it) {
    FullOperator op = base;
    const Real res = add_radiation(model, temp, lambda, op);
    reduce_operator(op, sys);
    const std::vector<Real> next = expand_solution(sys, solve_reduced(sys, kind));
    Real dmax = 0.0;
    for (std::size_t i = 0; i < next.size(); ++i)
      dmax = std::max(dmax, std::abs(next[i] - temp[i]));
    temp = next;
    if (!has_rad || (res < 1e-9 && dmax < 1e-9)) break;
  }
}

ThermalFields finalize(const Model& model, std::vector<Real> temperature) {
  ThermalFields res;
  res.flux_reaction = flux_reactions(model, temperature);
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
// (C/h) T_n added to the rhs, at step fraction `lambda` and step size `h`.
FullOperator step_operator(const Model& model, Index N, Real lambda, Real h,
                           const std::vector<Real>& temp) {
  const fem::FullThermalSystem fts = fem::assemble_full_thermal(model, lambda);
  FullOperator base = make_full_operator(N);
  seed_conduction(fts, base);
  for (std::size_t i = 0; i < fts.c_vals.size(); ++i) {
    const Real cij = fts.c_vals[i] / h;
    base.add(fts.c_rows[i], fts.c_cols[i], cij);
    base.b[static_cast<std::size_t>(fts.c_rows[i])] +=
        cij * temp[static_cast<std::size_t>(fts.c_cols[i])];
  }
  return base;
}

// Transient backward-Euler: (C/dt + Kt_film) T_{n+1} = (C/dt) T_n + q(+radiation),
// stepped from the initial field over the step period (increment.total) with the
// step's initial increment. Radiation is re-linearized inside each step.
ThermalFields solve_transient(const Model& model, LinearSystem& sys,
                              SolverKind kind) {
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
    const FullOperator base = step_operator(model, N, lambda, h, temp);
    nonlinear_solve(model, sys, base, kind, lambda, temp);
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

  if (model.procedure == Procedure::HeatTransferTransient)
    return solve_transient(model, sys, kind);

  // Steady state: conduction + film assembled once (full magnitude), radiation
  // iterated. Reuse the shared full assembly at lambda = 1; start from a zero field.
  const fem::FullThermalSystem fts = fem::assemble_full_thermal(model, 1.0);
  FullOperator base = make_full_operator(static_cast<Index>(model.mesh.num_nodes()));
  seed_conduction(fts, base);
  std::vector<Real> temp = expand_solution(
      sys, std::vector<Real>(static_cast<std::size_t>(sys.n_free), 0.0));
  nonlinear_solve(model, sys, base, kind, 1.0, temp);
  return finalize(model, std::move(temp));
}

CoupledFields solve_coupled(const Model& model, std::optional<SolverKind> forced) {
  // One-way (sequential) coupling. Step 1: solve the steady thermal field. The
  // coupled procedure carries the thermal cards (temperature BCs / flux / film), so
  // solve_heat_transfer's steady path applies verbatim (Coupled != transient).
  CoupledFields out;
  out.thermal = solve_heat_transfer(model, forced);

  // Step 2: copy the solved nodal temperatures into a mechanical copy of the model's
  // applied_temperature field (keyed by node id) and run the mechanical solve, which
  // adds the thermal-strain load and reports sigma = D (eps_mech - eps_th).
  Model mech = model;
  mech.procedure = Procedure::Static;
  mech.applied_temperature.clear();
  const Mesh& mesh = model.mesh;
  for (std::size_t ni = 0; ni < mesh.num_nodes(); ++ni)
    mech.applied_temperature[mesh.nodes()[ni].id] = out.thermal.temperature[ni];

  out.mechanical = solve_linear_static(mech, forced);
  return out;
}

}  // namespace cxpp::numerics
