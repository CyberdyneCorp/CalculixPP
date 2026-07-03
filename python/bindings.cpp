// Python bindings for CalculiX++ (spec: python-bindings).
// Exposes deck parsing and the linear-static solve, returning results as NumPy
// arrays for scripting and regression testing. C++ exceptions (incl. ParseError)
// propagate to Python as RuntimeError via pybind11's std::exception translation.
#include <optional>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <string>
#include <vector>

#include <map>

#include "calculixpp/compute/backend.hpp"
#include "calculixpp/core/element.hpp"
#include "calculixpp/io/inp_parser.hpp"
#include "calculixpp/io/results_writer.hpp"
#include "calculixpp/numerics/heat_transfer.hpp"
#include "calculixpp/numerics/linear_static.hpp"
#include "calculixpp/numerics/nonlinear_static.hpp"

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
  auto ri = ids.mutable_unchecked<1>();
  auto rc = coords.mutable_unchecked<2>();
  auto rt = temp.mutable_unchecked<1>();
  auto rr = rfl.mutable_unchecked<1>();
  for (std::size_t i = 0; i < n; ++i) {
    ri(static_cast<py::ssize_t>(i)) = m.mesh.nodes()[i].id;
    for (py::ssize_t j = 0; j < 3; ++j)
      rc(static_cast<py::ssize_t>(i), j) = m.mesh.nodes()[i].x[static_cast<std::size_t>(j)];
    rt(static_cast<py::ssize_t>(i)) = t.temperature[i];
    rr(static_cast<py::ssize_t>(i)) = t.flux_reaction[i];
  }
  py::dict d;
  d["node_ids"] = ids;
  d["node_coords"] = coords;
  d["temperature"] = temp;
  d["flux_reaction"] = rfl;
  d["num_nodes"] = n;
  d["num_elements"] = m.mesh.num_elements();
  return d;
}

py::dict solve_model(const Model& m, const std::string& solver,
                     const std::string& backend) {
  // Validate/select the backend first (may raise on an unknown name) and record
  // the one that actually ran. Phase 1 routes every solve through the CPU path.
  const std::string used = resolve_backend(backend);
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
  if (m.has_nonlinear_material()) {
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
