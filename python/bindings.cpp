// Python bindings for CalculiX++ (spec: python-bindings).
// Exposes deck parsing and the linear-static solve, returning results as NumPy
// arrays for scripting and regression testing. C++ exceptions (incl. ParseError)
// propagate to Python as RuntimeError via pybind11's std::exception translation.
#include <optional>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <string>
#include <vector>

#include <map>

#include "calculixpp/compute/backend.hpp"
#include "calculixpp/core/element.hpp"
#include "calculixpp/fem/assembly.hpp"
#include "calculixpp/io/inp_parser.hpp"
#include "calculixpp/io/results_writer.hpp"
#include "calculixpp/numerics/buckling.hpp"
#include "calculixpp/numerics/eigensolution.hpp"
#include "calculixpp/numerics/heat_transfer.hpp"
#include "calculixpp/numerics/direct_dynamics.hpp"
#include "calculixpp/numerics/modal_dynamics.hpp"
#include "calculixpp/numerics/linear_static.hpp"
#include "calculixpp/numerics/multistep.hpp"
#include "calculixpp/numerics/nonlinear_static.hpp"
#include "calculixpp/numerics/substructure.hpp"

namespace py = pybind11;
using namespace cxpp;

namespace {

template <std::size_t N>
py::array_t<double> rows_to_array(const std::vector<std::array<Real, N>>& v) {
  py::array_t<double> a({static_cast<py::ssize_t>(v.size()), static_cast<py::ssize_t>(N)});
  auto r = a.mutable_unchecked<2>();
  for (py::ssize_t i = 0; i < static_cast<py::ssize_t>(v.size()); ++i)
    for (py::ssize_t j = 0; j < static_cast<py::ssize_t>(N); ++j)
      r(i, j) = v[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)];
  return a;
}

py::array_t<double> vec3_to_array(const std::vector<Vec3>& v) {
  return rows_to_array<3>(v);
}

// Resolve the requested solver: an empty string honors the deck's SOLVER= (Auto
// by default, size-based); "cg"/"direct" force the path explicitly.
std::optional<numerics::SolverKind> kind_of(const std::string& s) {
  if (s.empty()) return std::nullopt;
  return s == "cg" ? numerics::SolverKind::CG : numerics::SolverKind::Direct;
}

std::string requested_solver_name(RequestedSolver r) {
  switch (r) {
    case RequestedSolver::Direct: return "direct";
    case RequestedSolver::CG: return "cg";
    default: return "auto";
  }
}

// Resolve and record the compute backend actually used for a solve. An empty
// string keeps the default (CPU); any known name is parsed and run through
// select_backend(), which in Phase 1 falls back to CPU for unimplemented targets.
// Returns the name of the backend that will actually run the solve. Throws
// (std::invalid_argument -> Python ValueError) on an unknown backend name.
std::string resolve_backend(const std::string& requested) {
  const compute::BackendKind want =
      requested.empty() ? compute::BackendKind::CPU : compute::backend_kind(requested);
  return compute::backend_name(compute::select_backend(want).kind());
}

// List the compute backends implemented on this build (CPU only in Phase 1).
std::vector<std::string> available_backend_names() {
  std::vector<std::string> names;
  for (const compute::BackendKind k : compute::available_backends())
    names.push_back(compute::backend_name(k));
  return names;
}

// Count solid elements by TYPE= name (C3D8, C3D20R, ...) for the summary. Only
// isoparametric mesh elements are counted; connectors are reported separately.
py::dict element_type_counts(const Model& m) {
  std::map<std::string, std::size_t> counts;  // ordered for a stable dict
  for (const Element& e : m.mesh.elements())
    counts[std::string(element_type_name(e.type))] += 1;
  py::dict d;
  for (const auto& [name, n] : counts) d[py::str(name)] = n;
  return d;
}

// Total number of higher-level constraint cards on the model (raw *EQUATION plus
// each *MPC / *RIGID BODY / *COUPLING / *TIE). This is the number of constraint
// CARDS, not the expanded linear-equation count (expand_constraints() would give
// the latter but resolves sets against the mesh and can throw).
std::size_t num_constraints(const Model& m) {
  return m.equations.size() + m.mpcs.size() + m.rigid_bodies.size() +
         m.couplings.size() + m.ties.size();
}

// A lightweight, solve-free view of a parsed deck. Beyond the Phase-1 counts it
// reflects the Phase-2 capabilities present in the deck (element-type breakdown,
// nonlinear-material / constraint / connector flags) so a script can inspect what
// a deck exercises — and which solve path it will take — without running a solve.
py::dict summary_dict(const Model& m) {
  std::vector<std::string> materials;
  materials.reserve(m.materials.size());
  for (const auto& [name, mat] : m.materials) materials.push_back(name);
  std::sort(materials.begin(), materials.end());
  py::dict d;
  d["num_nodes"] = m.mesh.num_nodes();
  d["num_elements"] = m.mesh.num_elements();
  d["num_materials"] = m.materials.size();
  d["materials"] = materials;
  d["requested_solver"] = requested_solver_name(m.solver);
  // Analysis procedure: "static" (mechanical) or heat transfer (steady/transient).
  d["procedure"] =
      m.procedure == Procedure::HeatTransferSteady      ? "heat transfer steady state"
      : m.procedure == Procedure::HeatTransferTransient ? "heat transfer transient"
      : m.procedure == Procedure::Coupled ? "coupled temperature-displacement"
      : m.procedure == Procedure::Frequency ? "frequency"
      : m.procedure == Procedure::Buckling ? "buckling"
      : m.procedure == Procedure::ModalDynamic ? "modal dynamic"
      : m.procedure == Procedure::SteadyStateDynamics ? "steady state dynamics"
      : m.procedure == Procedure::ComplexFrequency ? "complex frequency"
      : m.procedure == Procedure::Dynamic ? "dynamic"
      : m.procedure == Procedure::Substructure ? "substructure generate"
                                                        : "static";

  // Phase-2 introspection.
  d["element_type_counts"] = element_type_counts(m);
  d["has_plasticity"] = m.has_plasticity();
  // True when the deck routes to the incremental Newton driver (any *PLASTIC /
  // *HYPERELASTIC / *USER MATERIAL) — i.e. solve() will auto-dispatch to it.
  d["has_nonlinear_material"] = m.has_nonlinear_material();
  d["num_constraints"] = num_constraints(m);
  d["has_constraints"] = num_constraints(m) != 0;
  d["num_equations"] = m.equations.size();
  d["num_mpcs"] = m.mpcs.size();
  d["num_rigid_bodies"] = m.rigid_bodies.size();
  d["num_couplings"] = m.couplings.size();
  d["num_ties"] = m.ties.size();
  // Discrete/connector element counts (springs contribute to the static solve;
  // masses/dashpots are stored for dynamics).
  d["num_springs"] = m.springs.size();
  d["num_point_masses"] = m.point_masses.size();
  d["num_dashpots"] = m.dashpots.size();
  d["num_amplitudes"] = m.amplitudes.size();
  d["num_body_loads"] = m.body_loads.size();
  return d;
}

py::dict result_dict(const Model& m, const StaticFields& f) {
  const std::size_t n = m.mesh.num_nodes();
  py::array_t<int> ids(static_cast<py::ssize_t>(n));
  py::array_t<double> coords({static_cast<py::ssize_t>(n), py::ssize_t{3}});
  auto ri = ids.mutable_unchecked<1>();
  auto rc = coords.mutable_unchecked<2>();
  for (std::size_t i = 0; i < n; ++i) {
    ri(static_cast<py::ssize_t>(i)) = m.mesh.nodes()[i].id;
    for (py::ssize_t j = 0; j < 3; ++j)
      rc(static_cast<py::ssize_t>(i), j) = m.mesh.nodes()[i].x[static_cast<std::size_t>(j)];
  }
  py::dict d;
  d["node_ids"] = ids;
  d["node_coords"] = coords;
  d["displacement"] = vec3_to_array(f.displacement);
  d["stress"] = rows_to_array<6>(f.stress);
  d["strain"] = rows_to_array<6>(f.strain);
  d["reaction"] = vec3_to_array(f.reaction);
  d["num_nodes"] = n;
  d["num_elements"] = m.mesh.num_elements();
  // Contact results (CSTR) per slave node — a list of dicts (node, closed, pressure,
  // gap, tau). Empty for a contact-free deck. (spec: contact — contact output.)
  py::list contact;
  for (const ContactPoint& c : f.contact) {
    py::dict cd;
    cd["node"] = c.node_id;
    cd["closed"] = c.closed;
    cd["pressure"] = c.p;
    cd["gap"] = c.gap;
    cd["tau"] = c.tau;
    contact.append(cd);
  }
  d["contact"] = contact;
  return d;
}

// A heat-transfer result as NumPy arrays: node ids/coords plus temperature (NT) and
// heat-flux reaction (RFL), one scalar per node. (spec: heat-transfer-analysis.)
py::dict thermal_result_dict(const Model& m, const ThermalFields& t) {
  const std::size_t n = m.mesh.num_nodes();
  py::array_t<int> ids(static_cast<py::ssize_t>(n));
  py::array_t<double> coords({static_cast<py::ssize_t>(n), py::ssize_t{3}});
  py::array_t<double> temp(static_cast<py::ssize_t>(n));
  py::array_t<double> rfl(static_cast<py::ssize_t>(n));
  py::array_t<double> hfl({static_cast<py::ssize_t>(n), py::ssize_t{3}});  // nodal HFL
  auto ri = ids.mutable_unchecked<1>();
  auto rc = coords.mutable_unchecked<2>();
  auto rt = temp.mutable_unchecked<1>();
  auto rr = rfl.mutable_unchecked<1>();
  auto rh = hfl.mutable_unchecked<2>();
  for (std::size_t i = 0; i < n; ++i) {
    ri(static_cast<py::ssize_t>(i)) = m.mesh.nodes()[i].id;
    for (py::ssize_t j = 0; j < 3; ++j)
      rc(static_cast<py::ssize_t>(i), j) = m.mesh.nodes()[i].x[static_cast<std::size_t>(j)];
    rt(static_cast<py::ssize_t>(i)) = t.temperature[i];
    rr(static_cast<py::ssize_t>(i)) = t.flux_reaction[i];
    for (py::ssize_t j = 0; j < 3; ++j)
      rh(static_cast<py::ssize_t>(i), j) =
          t.heat_flux.empty() ? 0.0 : t.heat_flux[i][static_cast<std::size_t>(j)];
  }
  // Integration-point HFL (*EL PRINT HFL): elem ids, gp indices, and the flux rows.
  const std::size_t np = t.hfl_points.size();
  py::array_t<int> hp_elem(static_cast<py::ssize_t>(np));
  py::array_t<int> hp_gp(static_cast<py::ssize_t>(np));
  py::array_t<double> hp_flux({static_cast<py::ssize_t>(np), py::ssize_t{3}});
  auto pe = hp_elem.mutable_unchecked<1>();
  auto pg = hp_gp.mutable_unchecked<1>();
  auto pf = hp_flux.mutable_unchecked<2>();
  for (std::size_t i = 0; i < np; ++i) {
    pe(static_cast<py::ssize_t>(i)) = t.hfl_points[i].elem_id;
    pg(static_cast<py::ssize_t>(i)) = t.hfl_points[i].gp;
    for (py::ssize_t j = 0; j < 3; ++j)
      pf(static_cast<py::ssize_t>(i), j) = t.hfl_points[i].flux[static_cast<std::size_t>(j)];
  }
  py::dict d;
  d["node_ids"] = ids;
  d["node_coords"] = coords;
  d["temperature"] = temp;
  d["flux_reaction"] = rfl;
  d["heat_flux"] = hfl;
  d["hfl_elem"] = hp_elem;
  d["hfl_gp"] = hp_gp;
  d["hfl_flux"] = hp_flux;
  d["num_nodes"] = n;
  d["num_elements"] = m.mesh.num_elements();
  return d;
}

// A *FREQUENCY result as NumPy arrays: eigenvalues (λ = ω²), natural angular
// frequencies (ω) and cyclic frequencies (f = ω/2π), the mass-normalized mode shapes
// (n_modes x n_nodes x 3), and the modal participation factors / effective mass for
// the three translational excitation directions. (spec: eigensolution, modal-and-
// buckling — *FREQUENCY.)
py::dict frequency_result_dict(const Model& m, const numerics::EigenBasis& basis,
                               const fem::LinearSystem& M) {
  const std::size_t nm = basis.modes.size();
  const std::size_t nn = m.mesh.num_nodes();
  py::array_t<double> eigv(static_cast<py::ssize_t>(nm));
  py::array_t<double> omega(static_cast<py::ssize_t>(nm));
  py::array_t<double> freq(static_cast<py::ssize_t>(nm));
  py::array_t<double> shapes({static_cast<py::ssize_t>(nm),
                              static_cast<py::ssize_t>(nn), py::ssize_t{3}});
  auto re = eigv.mutable_unchecked<1>();
  auto ro = omega.mutable_unchecked<1>();
  auto rf = freq.mutable_unchecked<1>();
  auto rs = shapes.mutable_unchecked<3>();
  for (std::size_t k = 0; k < nm; ++k) {
    re(static_cast<py::ssize_t>(k)) = basis.modes[k].eigenvalue;
    ro(static_cast<py::ssize_t>(k)) = basis.modes[k].omega;
    rf(static_cast<py::ssize_t>(k)) = basis.modes[k].frequency;
    for (std::size_t i = 0; i < nn; ++i)
      for (py::ssize_t j = 0; j < 3; ++j)
        rs(static_cast<py::ssize_t>(k), static_cast<py::ssize_t>(i), j) =
            basis.modes[k].shape[i][static_cast<std::size_t>(j)];
  }
  // Participation factors / effective mass (n_modes x 3) and totals per direction.
  py::array_t<double> part({static_cast<py::ssize_t>(nm), py::ssize_t{3}});
  py::array_t<double> effmass({static_cast<py::ssize_t>(nm), py::ssize_t{3}});
  py::array_t<double> total_eff(py::ssize_t{3});
  auto rp = part.mutable_unchecked<2>();
  auto rem = effmass.mutable_unchecked<2>();
  auto rte = total_eff.mutable_unchecked<1>();
  for (int dir = 0; dir < 3; ++dir) {
    const numerics::Participation p = numerics::participation(basis, M, dir);
    for (std::size_t k = 0; k < nm; ++k) {
      rp(static_cast<py::ssize_t>(k), dir) = p.factor[k];
      rem(static_cast<py::ssize_t>(k), dir) = p.effective_mass[k];
    }
    rte(dir) = p.total_effective_mass;
  }
  py::dict d;
  d["procedure"] = "frequency";
  d["num_modes"] = nm;
  d["num_nodes"] = nn;
  d["num_elements"] = m.mesh.num_elements();
  d["eigenvalue"] = eigv;
  d["omega"] = omega;
  d["frequency"] = freq;
  d["mode_shape"] = shapes;
  d["participation"] = part;
  d["effective_mass"] = effmass;
  d["total_effective_mass"] = total_eff;
  return d;
}

// A *COMPLEX FREQUENCY deck returns the damped complex modes: complex eigenvalues (split
// into real/imag arrays), damped frequencies f_d, damping ratios ζ, undamped angular
// frequencies ω_n, and the real+imag parts of each complex mode shape (n_modes x n_nodes
// x 3). Option-(B) proportional damping (spec: modal-and-buckling — complex frequency).
py::dict complex_frequency_result_dict(const Model& m,
                                       const numerics::ComplexEigenBasis& cx) {
  const std::size_t nm = cx.modes.size();
  const std::size_t nn = m.mesh.num_nodes();
  py::array_t<double> eig_re(static_cast<py::ssize_t>(nm));
  py::array_t<double> eig_im(static_cast<py::ssize_t>(nm));
  py::array_t<double> fd(static_cast<py::ssize_t>(nm));
  py::array_t<double> zeta(static_cast<py::ssize_t>(nm));
  py::array_t<double> omega_n(static_cast<py::ssize_t>(nm));
  py::array_t<double> shape_re({static_cast<py::ssize_t>(nm),
                                static_cast<py::ssize_t>(nn), py::ssize_t{3}});
  py::array_t<double> shape_im({static_cast<py::ssize_t>(nm),
                                static_cast<py::ssize_t>(nn), py::ssize_t{3}});
  auto rer = eig_re.mutable_unchecked<1>();
  auto rei = eig_im.mutable_unchecked<1>();
  auto rfd = fd.mutable_unchecked<1>();
  auto rz = zeta.mutable_unchecked<1>();
  auto ron = omega_n.mutable_unchecked<1>();
  auto rsr = shape_re.mutable_unchecked<3>();
  auto rsi = shape_im.mutable_unchecked<3>();
  for (std::size_t k = 0; k < nm; ++k) {
    const numerics::ComplexMode& md = cx.modes[k];
    const auto kk = static_cast<py::ssize_t>(k);
    rer(kk) = md.eigenvalue.real();
    rei(kk) = md.eigenvalue.imag();
    rfd(kk) = md.frequency;
    rz(kk) = md.zeta;
    ron(kk) = md.omega_n;
    for (std::size_t i = 0; i < nn; ++i)
      for (py::ssize_t j = 0; j < 3; ++j) {
        rsr(kk, static_cast<py::ssize_t>(i), j) = md.shape_real[i][static_cast<std::size_t>(j)];
        rsi(kk, static_cast<py::ssize_t>(i), j) = md.shape_imag[i][static_cast<std::size_t>(j)];
      }
  }
  py::dict d;
  d["procedure"] = "complex frequency";
  d["num_modes"] = nm;
  d["num_nodes"] = nn;
  d["num_elements"] = m.mesh.num_elements();
  d["eigenvalues_real"] = eig_re;
  d["eigenvalues_imag"] = eig_im;
  d["damped_frequencies"] = fd;
  d["damping_ratios"] = zeta;
  d["omega_n"] = omega_n;
  d["mode_shapes_real"] = shape_re;
  d["mode_shapes_imag"] = shape_im;
  return d;
}

// A *BUCKLE result as NumPy arrays: the buckling load factors λ (ascending positive)
// and the buckling mode shapes (n_modes x n_nodes x 3). The critical load is
// factors[0] * f_ref. (spec: modal-and-buckling — *BUCKLE.)
py::dict buckling_result_dict(const Model& m, const numerics::BucklingReport& rep) {
  const std::size_t nm = rep.basis.modes.size();
  const std::size_t nn = m.mesh.num_nodes();
  py::array_t<double> factors(static_cast<py::ssize_t>(nm));
  py::array_t<double> shapes({static_cast<py::ssize_t>(nm),
                              static_cast<py::ssize_t>(nn), py::ssize_t{3}});
  auto rfac = factors.mutable_unchecked<1>();
  auto rs = shapes.mutable_unchecked<3>();
  for (std::size_t k = 0; k < nm; ++k) {
    rfac(static_cast<py::ssize_t>(k)) = rep.basis.modes[k].eigenvalue;
    for (std::size_t i = 0; i < nn; ++i)
      for (py::ssize_t j = 0; j < 3; ++j)
        rs(static_cast<py::ssize_t>(k), static_cast<py::ssize_t>(i), j) =
            rep.basis.modes[k].shape[i][static_cast<std::size_t>(j)];
  }
  py::dict d;
  d["procedure"] = "buckling";
  d["num_modes"] = nm;
  d["num_nodes"] = nn;
  d["num_elements"] = m.mesh.num_elements();
  d["factors"] = factors;
  d["mode_shape"] = shapes;
  return d;
}

// A *SUBSTRUCTURE GENERATE deck condenses onto the retained DOFs (Craig-Bampton /
// Guyan) and returns the reduced operators + their labels. `k_reduced`/`m_reduced` are
// dim×dim arrays (dim = retained DOFs + fixed-interface modes); `retained_node` /
// `retained_comp` name the leading retained block; `modal_omega` the appended modal
// block's frequencies. (spec: substructure-generation — reachable from Python.)
py::dict substructure_result_dict(const Model& m) {
  const numerics::Superelement se = numerics::generate_substructure(m);
  const std::size_t dim = se.dim();
  py::array_t<double> K({static_cast<py::ssize_t>(dim),
                         static_cast<py::ssize_t>(dim)});
  auto rk = K.mutable_unchecked<2>();
  for (std::size_t i = 0; i < dim; ++i)
    for (std::size_t j = 0; j < dim; ++j)
      rk(static_cast<py::ssize_t>(i), static_cast<py::ssize_t>(j)) =
          se.k_reduced[i * dim + j];
  py::dict d;
  d["procedure"] = "substructure generate";
  d["num_retained"] = se.n_retained;
  d["num_modes"] = se.n_modes;
  d["dim"] = dim;
  d["k_reduced"] = K;
  if (!se.m_reduced.empty()) {
    py::array_t<double> M({static_cast<py::ssize_t>(dim),
                           static_cast<py::ssize_t>(dim)});
    auto rm = M.mutable_unchecked<2>();
    for (std::size_t i = 0; i < dim; ++i)
      for (std::size_t j = 0; j < dim; ++j)
        rm(static_cast<py::ssize_t>(i), static_cast<py::ssize_t>(j)) =
            se.m_reduced[i * dim + j];
    d["m_reduced"] = M;
  }
  d["retained_node"] = se.retained_node;
  d["retained_comp"] = se.retained_comp;
  d["modal_omega"] = se.modal_omega;
  return d;
}

// Build the modal system (eigenbasis + damping projection) shared by *MODAL DYNAMIC and
// *STEADY STATE DYNAMICS. Extracts the requested number of modes, applies the deck's
// Rayleigh / modal damping, and returns the reduced modal operators plus the free-DOF
// load pattern (K.rhs — the reduced concentrated-load vector) for superposition. The
// eigenbasis is returned by value into `basis_out` because ModalSystem borrows it.
numerics::ModalSystem build_modal_system(const Model& m,
                                         numerics::EigenBasis& basis_out,
                                         std::vector<Real>& load_out,
                                         const fem::LinearSystem& K,
                                         const fem::LinearSystem& M) {
  const std::size_t nreq = m.num_eigenvalues > 0
                               ? static_cast<std::size_t>(m.num_eigenvalues)
                               : static_cast<std::size_t>(K.n_free);
  basis_out = numerics::extract_modes(K, M, nreq);
  numerics::Damping damp;
  damp.alpha = m.rayleigh.alpha;
  damp.beta = m.rayleigh.beta;
  damp.modal_ratios = m.modal_damping;
  load_out = K.rhs;  // reduced free-DOF concentrated-load pattern
  return numerics::project_modal_system(basis_out, damp);
}

// A *MODAL DYNAMIC result: the transient nodal displacement history (n_steps x n_nodes
// x 3) sampled at the step times, plus the modal frequencies/damping used. (spec:
// dynamic-analysis — modal dynamic.)
py::dict modal_dynamic_result_dict(const Model& m) {
  const fem::LinearSystem K = fem::assemble_linear_static(m);
  const fem::LinearSystem M = fem::assemble_mass(m, /*lumped=*/false);
  numerics::EigenBasis basis;
  std::vector<Real> pattern;
  const numerics::ModalSystem sys = build_modal_system(m, basis, pattern, K, M);

  numerics::ModalLoad load;
  load.pattern = pattern;
  const Real dt = m.dynamic_dt > 0.0 ? m.dynamic_dt : 1.0;
  const Real t_end = m.dynamic_t_end > 0.0 ? m.dynamic_t_end : dt;
  const std::vector<numerics::ModalTimePoint> hist =
      numerics::modal_dynamic(sys, load, dt, t_end);

  const std::size_t nt = hist.size(), nn = m.mesh.num_nodes();
  py::array_t<double> times(static_cast<py::ssize_t>(nt));
  py::array_t<double> disp({static_cast<py::ssize_t>(nt),
                            static_cast<py::ssize_t>(nn), py::ssize_t{3}});
  auto rt = times.mutable_unchecked<1>();
  auto rd = disp.mutable_unchecked<3>();
  for (std::size_t s = 0; s < nt; ++s) {
    rt(static_cast<py::ssize_t>(s)) = hist[s].time;
    for (std::size_t i = 0; i < nn; ++i)
      for (py::ssize_t j = 0; j < 3; ++j)
        rd(static_cast<py::ssize_t>(s), static_cast<py::ssize_t>(i), j) =
            hist[s].displacement[i][static_cast<std::size_t>(j)];
  }
  py::array_t<double> zeta(static_cast<py::ssize_t>(sys.n_modes));
  py::array_t<double> omega(static_cast<py::ssize_t>(sys.n_modes));
  auto rz = zeta.mutable_unchecked<1>();
  auto ro = omega.mutable_unchecked<1>();
  for (std::size_t k = 0; k < sys.n_modes; ++k) {
    rz(static_cast<py::ssize_t>(k)) = sys.zeta[k];
    ro(static_cast<py::ssize_t>(k)) = sys.omega[k];
  }
  py::dict d;
  d["procedure"] = "modal dynamic";
  d["num_modes"] = sys.n_modes;
  d["num_nodes"] = nn;
  d["num_steps"] = nt;
  d["time"] = times;
  d["displacement"] = disp;  // (n_steps, n_nodes, 3)
  d["omega"] = omega;
  d["zeta"] = zeta;
  return d;
}

// A *STEADY STATE DYNAMICS result: the harmonic response amplitude and phase per node
// over the frequency sweep (n_freq x n_nodes x 3 each), plus the swept frequencies.
// (spec: dynamic-analysis — steady-state dynamics.)
py::dict steady_state_result_dict(const Model& m) {
  const fem::LinearSystem K = fem::assemble_linear_static(m);
  const fem::LinearSystem M = fem::assemble_mass(m, /*lumped=*/false);
  numerics::EigenBasis basis;
  std::vector<Real> pattern;
  const numerics::ModalSystem sys = build_modal_system(m, basis, pattern, K, M);

  const Real f_lo = m.steady_f_lo;
  const Real f_hi = m.steady_f_hi > 0.0 ? m.steady_f_hi : f_lo;
  const std::size_t npts =
      m.steady_num_points > 0 ? static_cast<std::size_t>(m.steady_num_points) : 20;
  const std::vector<numerics::HarmonicResponse> sweep =
      numerics::steady_state_sweep(sys, pattern, f_lo, f_hi, npts);

  const std::size_t nf = sweep.size(), nn = m.mesh.num_nodes();
  py::array_t<double> freqs(static_cast<py::ssize_t>(nf));
  py::array_t<double> ampl({static_cast<py::ssize_t>(nf),
                            static_cast<py::ssize_t>(nn), py::ssize_t{3}});
  py::array_t<double> phase({static_cast<py::ssize_t>(nf),
                            static_cast<py::ssize_t>(nn), py::ssize_t{3}});
  auto rff = freqs.mutable_unchecked<1>();
  auto ra = ampl.mutable_unchecked<3>();
  auto rp = phase.mutable_unchecked<3>();
  // The complex free-DOF amplitude is expanded to full nodal components through the
  // constraint transform (SPC value 0 — a harmonic response has no prescribed
  // displacement), separately for the real and imaginary parts, then converted to
  // magnitude/phase per node component.
  const std::vector<Real> zero_prescribed(K.prescribed.size(), 0.0);
  for (std::size_t s = 0; s < nf; ++s) {
    rff(static_cast<py::ssize_t>(s)) = sweep[s].frequency;
    std::vector<Real> re(sweep[s].amplitude.size()), im(sweep[s].amplitude.size());
    for (std::size_t i = 0; i < re.size(); ++i) {
      re[i] = sweep[s].amplitude[i].real();
      im[i] = sweep[s].amplitude[i].imag();
    }
    for (std::size_t i = 0; i < nn; ++i)
      for (int c = 0; c < 3; ++c) {
        const std::size_t g = i * static_cast<std::size_t>(kDofsPerNode) +
                              static_cast<std::size_t>(c);
        const Real ur = K.transform.displacement(g, re, zero_prescribed);
        const Real ui = K.transform.displacement(g, im, zero_prescribed);
        ra(static_cast<py::ssize_t>(s), static_cast<py::ssize_t>(i), c) =
            std::sqrt(ur * ur + ui * ui);
        rp(static_cast<py::ssize_t>(s), static_cast<py::ssize_t>(i), c) =
            std::atan2(ui, ur);
      }
  }
  py::dict d;
  d["procedure"] = "steady state dynamics";
  d["num_modes"] = sys.n_modes;
  d["num_nodes"] = nn;
  d["num_frequencies"] = nf;
  d["frequency"] = freqs;
  d["amplitude"] = ampl;  // (n_freq, n_nodes, 3)
  d["phase"] = phase;     // (n_freq, n_nodes, 3) radians
  return d;
}

// A *DYNAMIC result: the direct HHT-α time integration of M a + C v + K u = f(t). Returns
// the full displacement/velocity/acceleration history (n_steps x n_nodes x 3) plus the
// per-step energy (kinetic/strain/total) and the energy drift over the run. (spec:
// dynamic-analysis — direct-integration dynamics.)
py::dict dynamic_result_dict(const Model& m) {
  numerics::DirectOptions opts;
  opts.hht = numerics::HhtParams::from_alpha(m.dynamic_alpha);
  opts.damping.alpha = m.rayleigh.alpha;
  opts.damping.beta = m.rayleigh.beta;
  opts.nonlinear = m.dynamic_nonlinear;
  const Real dt = m.dynamic_dt > 0.0 ? m.dynamic_dt : 1.0;
  const Real t_end = m.dynamic_t_end > 0.0 ? m.dynamic_t_end : dt;
  numerics::DirectReport rep;
  const std::vector<numerics::DirectTimePoint> hist =
      numerics::direct_dynamic(m, dt, t_end, opts, nullptr, &rep);

  const std::size_t nt = hist.size(), nn = m.mesh.num_nodes();
  py::array_t<double> times(static_cast<py::ssize_t>(nt));
  py::array_t<double> disp({static_cast<py::ssize_t>(nt),
                            static_cast<py::ssize_t>(nn), py::ssize_t{3}});
  py::array_t<double> vel({static_cast<py::ssize_t>(nt),
                           static_cast<py::ssize_t>(nn), py::ssize_t{3}});
  py::array_t<double> acc({static_cast<py::ssize_t>(nt),
                           static_cast<py::ssize_t>(nn), py::ssize_t{3}});
  py::array_t<double> energy(static_cast<py::ssize_t>(nt));
  auto rt = times.mutable_unchecked<1>();
  auto rd = disp.mutable_unchecked<3>();
  auto rv = vel.mutable_unchecked<3>();
  auto ra = acc.mutable_unchecked<3>();
  auto re = energy.mutable_unchecked<1>();
  for (std::size_t s = 0; s < nt; ++s) {
    const auto ss = static_cast<py::ssize_t>(s);
    rt(ss) = hist[s].time;
    re(ss) = hist[s].total_energy;
    for (std::size_t i = 0; i < nn; ++i)
      for (py::ssize_t j = 0; j < 3; ++j) {
        const auto jj = static_cast<std::size_t>(j);
        const auto ii = static_cast<py::ssize_t>(i);
        rd(ss, ii, j) = hist[s].displacement[i][jj];
        rv(ss, ii, j) = hist[s].velocity[i][jj];
        ra(ss, ii, j) = hist[s].acceleration[i][jj];
      }
  }
  py::dict d;
  d["procedure"] = "dynamic";
  d["num_nodes"] = nn;
  d["num_steps"] = nt;
  d["hht_alpha"] = opts.hht.alpha;
  d["nonlinear"] = rep.nonlinear;
  d["newton_iterations"] = rep.iterations;
  d["energy_drift"] = rep.energy_drift;
  d["time"] = times;
  d["displacement"] = disp;    // (n_steps, n_nodes, 3)
  d["velocity"] = vel;         // (n_steps, n_nodes, 3)
  d["acceleration"] = acc;     // (n_steps, n_nodes, 3)
  d["total_energy"] = energy;  // (n_steps,)
  return d;
}

py::dict solve_model(const Model& m, const std::string& solver,
                     const std::string& backend) {
  // Validate/select the backend first (may raise on an unknown name) and record
  // the one that actually ran. Phase 1 routes every solve through the CPU path.
  const std::string used = resolve_backend(backend);
  // A *FREQUENCY deck extracts the lowest N natural frequencies / mode shapes of the
  // generalized eigenproblem K x = λ M x. (spec: modal-and-buckling — *FREQUENCY.)
  if (m.procedure == Procedure::Frequency) {
    const std::size_t nreq = m.num_eigenvalues > 0
                                 ? static_cast<std::size_t>(m.num_eigenvalues)
                                 : 1;
    const fem::LinearSystem K = fem::assemble_linear_static(m);
    const fem::LinearSystem M = fem::assemble_mass(m, /*lumped=*/false);
    const numerics::EigenBasis basis = numerics::extract_modes(K, M, nreq);
    py::dict d = frequency_result_dict(m, basis, M);
    d["backend"] = used;
    return d;
  }
  // A *BUCKLE deck runs the two-step prestress driver and returns the buckling factors
  // (ascending positive) + mode shapes. (spec: modal-and-buckling — *BUCKLE.)
  if (m.procedure == Procedure::Buckling) {
    const std::size_t nreq = m.num_buckling_modes > 0
                                 ? static_cast<std::size_t>(m.num_buckling_modes)
                                 : 1;
    const numerics::BucklingReport rep = numerics::solve_buckling(m, nreq);
    py::dict d = buckling_result_dict(m, rep);
    d["backend"] = used;
    return d;
  }
  // A *MODAL DYNAMIC deck integrates the decoupled modal SDOFs over the *FREQUENCY basis
  // and returns the transient displacement history. (spec: dynamic-analysis.)
  if (m.procedure == Procedure::ModalDynamic) {
    py::dict d = modal_dynamic_result_dict(m);
    d["backend"] = used;
    return d;
  }
  // A *STEADY STATE DYNAMICS deck sweeps the harmonic response over the requested
  // frequency band by modal superposition. (spec: dynamic-analysis.)
  if (m.procedure == Procedure::SteadyStateDynamics) {
    py::dict d = steady_state_result_dict(m);
    d["backend"] = used;
    return d;
  }
  // A *COMPLEX FREQUENCY deck reduces the proportional damping onto the *FREQUENCY basis
  // and returns damped complex modes (eigenvalues, damped frequencies, damping ratios,
  // complex mode shapes). (spec: modal-and-buckling — complex frequency; option B.)
  if (m.procedure == Procedure::ComplexFrequency) {
    const std::size_t nreq = m.num_complex_modes > 0
                                 ? static_cast<std::size_t>(m.num_complex_modes)
                                 : 1;
    const fem::LinearSystem K = fem::assemble_linear_static(m);
    const fem::LinearSystem M = fem::assemble_mass(m, /*lumped=*/false);
    const numerics::EigenBasis basis = numerics::extract_modes(K, M, nreq);
    numerics::Damping damp;
    damp.alpha = m.rayleigh.alpha;
    damp.beta = m.rayleigh.beta;
    damp.modal_ratios = m.modal_damping;
    const numerics::ComplexEigenBasis cx =
        numerics::extract_complex_modes(basis, damp, nreq);
    py::dict d = complex_frequency_result_dict(m, cx);
    d["backend"] = used;
    return d;
  }
  // A *DYNAMIC deck runs direct HHT-α time integration in physical coordinates and
  // returns the displacement/velocity/acceleration + energy history. (spec:
  // dynamic-analysis — direct-integration dynamics.)
  if (m.procedure == Procedure::Dynamic) {
    py::dict d = dynamic_result_dict(m);
    d["backend"] = used;
    return d;
  }
  // A *SUBSTRUCTURE GENERATE deck condenses the model onto the retained DOFs and returns
  // the reduced stiffness (+ mass, Craig-Bampton). (spec: substructure-generation.)
  if (m.procedure == Procedure::Substructure) {
    py::dict d = substructure_result_dict(m);
    d["backend"] = used;
    return d;
  }
  // A *COUPLED TEMPERATURE-DISPLACEMENT deck solves the thermal field, applies the
  // resulting thermal strain, then the mechanical field (one-way coupling). The
  // result dict carries the mechanical fields plus the temperature/flux_reaction.
  if (m.procedure == Procedure::Coupled) {
    const CoupledFields cf = numerics::solve_coupled(m, kind_of(solver));
    py::dict cd = result_dict(m, cf.mechanical);
    py::dict td = thermal_result_dict(m, cf.thermal);
    cd["temperature"] = td["temperature"];
    cd["flux_reaction"] = td["flux_reaction"];
    cd["backend"] = used;
    cd["procedure"] = "coupled temperature-displacement";
    return cd;
  }
  // A *HEAT TRANSFER deck solves the scalar temperature field; the mechanical
  // path is unchanged.
  if (m.procedure == Procedure::HeatTransferSteady ||
      m.procedure == Procedure::HeatTransferTransient) {
    py::dict td = thermal_result_dict(
        m, numerics::solve_heat_transfer(m, kind_of(solver)));
    td["backend"] = used;
    td["procedure"] = m.procedure == Procedure::HeatTransferTransient
                          ? "heat transfer transient"
                          : "heat transfer steady state";
    return td;
  }
  // A model with a nonlinear material (*PLASTIC / *HYPERELASTIC / *USER MATERIAL)
  // routes to the Newton-Raphson driver (load applied incrementally); a purely
  // linear-elastic model keeps the direct linear path unchanged.
  py::dict d;
  // A *CONTACT PAIR deck (nonlinear constraint) or a nonlinear material routes to the
  // Newton driver; a purely linear-elastic contact-free model keeps the linear path.
  if (m.has_nonlinear_material() || m.has_contact()) {
    numerics::NonlinearOptions opts;
    opts.forced = kind_of(solver);
    numerics::NonlinearReport rep;
    d = result_dict(m, numerics::solve_nonlinear_static(m, opts, &rep));
    d["newton_increments"] = rep.increments;
    d["newton_iterations"] = rep.iterations;
    d["newton_cutbacks"] = rep.cutbacks;
    d["converged"] = rep.converged;
  } else {
    d = result_dict(m, numerics::solve_linear_static(m, kind_of(solver)));
  }
  d["backend"] = used;
  return d;
}

py::dict solve_file(const std::string& path, const std::string& solver,
                    const std::string& backend) {
  return solve_model(io::parse_inp_file(path), solver, backend);
}

py::dict solve_text(const std::string& text, const std::string& solver,
                    const std::string& backend) {
  return solve_model(io::parse_inp(text), solver, backend);
}

// Nonlinear Newton-Raphson solve. Reproduces the linear solve on a linear model.
// The returned dict adds newton_increments / newton_iterations / converged so
// callers can inspect the solution path.
py::dict solve_nonlinear_model(const Model& m, const std::string& solver,
                               const std::string& backend, bool line_search) {
  const std::string used = resolve_backend(backend);
  numerics::NonlinearOptions opts;
  opts.line_search = line_search;
  opts.forced = kind_of(solver);
  numerics::NonlinearReport rep;
  py::dict d = result_dict(m, numerics::solve_nonlinear_static(m, opts, &rep));
  d["backend"] = used;
  d["newton_increments"] = rep.increments;
  d["newton_iterations"] = rep.iterations;
  d["newton_cutbacks"] = rep.cutbacks;
  d["converged"] = rep.converged;
  return d;
}

py::dict solve_nonlinear_file(const std::string& path, const std::string& solver,
                              const std::string& backend, bool line_search) {
  return solve_nonlinear_model(io::parse_inp_file(path), solver, backend,
                               line_search);
}

py::dict solve_nonlinear_text(const std::string& text, const std::string& solver,
                              const std::string& backend, bool line_search) {
  return solve_nonlinear_model(io::parse_inp(text), solver, backend, line_search);
}

// Multi-step linear-static solve (spec: multi-step analysis). Parses the deck into its
// per-*STEP models and runs the step-loop driver, returning a LIST of per-step result
// dicts (one per *STEP, each carrying the accumulated total displacement/stress/strain/
// reaction at the end of that step). A single-*STEP deck returns a one-element list
// whose fields equal solve()'s (the single-step fast path). Enables validating a stock
// two-step deck's per-step *NODE PRINT blocks. `solver`/`backend` as in solve().
py::list solve_multistep_models(const std::vector<Model>& steps,
                                const std::string& solver, const std::string& backend) {
  const std::string used = resolve_backend(backend);
  numerics::MultiStepReport rep;
  const std::vector<StaticFields> all =
      numerics::solve_multistep_static_all(steps, kind_of(solver), &rep);
  py::list out;
  for (std::size_t i = 0; i < all.size(); ++i) {
    py::dict d = result_dict(steps[i], all[i]);
    d["backend"] = used;
    d["step"] = static_cast<int>(i + 1);
    d["num_steps"] = static_cast<int>(all.size());
    out.append(d);
  }
  return out;
}

py::list solve_multistep_file(const std::string& path, const std::string& solver,
                              const std::string& backend) {
  return solve_multistep_models(io::parse_inp_steps_file(path), solver, backend);
}

py::list solve_multistep_text(const std::string& text, const std::string& solver,
                              const std::string& backend) {
  return solve_multistep_models(io::parse_inp_steps(text), solver, backend);
}

py::dict summary_file(const std::string& path) {
  return summary_dict(io::parse_inp_file(path));
}

py::dict summary_text(const std::string& text) {
  return summary_dict(io::parse_inp(text));
}

}  // namespace

PYBIND11_MODULE(calculixpp, mod) {
  mod.doc() = "CalculiX++ — modern C++20 finite element solver (Python bindings)";
  mod.def("solve", &solve_file, py::arg("path"), py::arg("solver") = "",
          py::arg("backend") = "",
          "Parse an .inp deck and solve the linear-static step. solver='' honors "
          "the deck's SOLVER= (default direct); pass 'direct' or 'cg' to force a "
          "path. backend='' uses the default compute backend; pass a name from "
          "available_backends() (unimplemented backends fall back to cpu). Returns "
          "a dict of NumPy arrays: node_ids, node_coords, displacement (Nx3), "
          "stress (Nx6, xx/yy/zz/xy/xz/yz), strain (Nx6, same order, engineering "
          "shear), reaction (Nx3), plus num_nodes, num_elements, and backend (the "
          "backend that actually ran the solve). A *HEAT TRANSFER, STEADY STATE deck "
          "auto-dispatches to the thermal solver and instead returns temperature "
          "(N,) and flux_reaction (N,) with procedure='heat transfer steady state'.");
  mod.def("solve_text", &solve_text, py::arg("text"), py::arg("solver") = "",
          py::arg("backend") = "",
          "Like solve() but takes the deck contents as a string.");

  mod.def("solve_multistep", &solve_multistep_file, py::arg("path"),
          py::arg("solver") = "", py::arg("backend") = "",
          "Parse a multi-*STEP .inp deck and solve every step in order with the "
          "step-loop driver (spec: multi-step analysis), carrying the converged state "
          "forward. Returns a LIST of per-step result dicts (same fields as solve(), "
          "plus 'step' 1-based and 'num_steps'); entry i is the accumulated total at "
          "the end of step i+1. Two linear load steps sum to the single-step result "
          "(superposition); *MODEL CHANGE removes/re-adds elements across steps "
          "(strain-free reactivation); *BOUNDARY,FIXED holds a DOF at its deformed "
          "value; OP=NEW/MOD accumulate loads and *CHANGE SOLID SECTION rebinds at the "
          "boundary. A single-*STEP deck returns a one-element list equal to solve().");
  mod.def("solve_multistep_text", &solve_multistep_text, py::arg("text"),
          py::arg("solver") = "", py::arg("backend") = "",
          "Like solve_multistep() but takes the deck contents as a string.");

  mod.def("solve_nonlinear", &solve_nonlinear_file, py::arg("path"),
          py::arg("solver") = "", py::arg("backend") = "",
          py::arg("line_search") = false,
          "Parse an .inp deck and solve the static step with the Newton-Raphson "
          "driver (spec: nonlinear-solution-control). On a linear model this "
          "reproduces solve() exactly. Returns the same fields as solve() plus "
          "newton_increments, newton_iterations, newton_cutbacks, and converged. "
          "line_search=True enables the optional Newton step scaling (off by "
          "default).");
  mod.def("solve_nonlinear_text", &solve_nonlinear_text, py::arg("text"),
          py::arg("solver") = "", py::arg("backend") = "",
          py::arg("line_search") = false,
          "Like solve_nonlinear() but takes the deck contents as a string.");

  mod.def("summary", &summary_file, py::arg("path"),
          "Parse an .inp deck without solving and return a dict describing it. "
          "Phase-1 fields: num_nodes, num_elements, num_materials, materials "
          "(sorted names), requested_solver ('auto'/'direct'/'cg'). Phase-2 "
          "introspection: element_type_counts (dict TYPE=->count), has_plasticity, "
          "has_nonlinear_material (True iff solve() auto-routes to the Newton "
          "driver), num_constraints/has_constraints and the per-kind breakdown "
          "(num_equations, num_mpcs, num_rigid_bodies, num_couplings, num_ties), "
          "connector counts (num_springs, num_point_masses, num_dashpots), "
          "num_amplitudes, and num_body_loads.");
  mod.def("summary_text", &summary_text, py::arg("text"),
          "Like summary() but takes the deck contents as a string.");

  mod.def("available_backends", &available_backend_names,
          "List the compute backends implemented on this build (names usable as "
          "the backend= argument to solve()). Phase 1 exposes only 'cpu'.");
  mod.def(
      "selected_backend",
      [](const std::string& requested) { return resolve_backend(requested); },
      py::arg("backend") = "",
      "Return the name of the compute backend that would run a solve for the "
      "given request. An empty request resolves the default ('cpu'); an "
      "unimplemented backend name resolves to its CPU fallback; an unknown name "
      "raises ValueError.");
}
