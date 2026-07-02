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

#include "calculixpp/compute/backend.hpp"
#include "calculixpp/io/inp_parser.hpp"
#include "calculixpp/io/results_writer.hpp"
#include "calculixpp/numerics/linear_static.hpp"

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

// A lightweight, solve-free view of a parsed deck: counts and material names.
// Useful for scripting/inspection without running (or being able to run) a solve.
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
  d["requested_solver"] =
      requested_solver_name(m.solver);
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

py::dict solve_model(const Model& m, const std::string& solver,
                     const std::string& backend) {
  // Validate/select the backend first (may raise on an unknown name) and record
  // the one that actually ran. Phase 1 routes every solve through the CPU path.
  const std::string used = resolve_backend(backend);
  py::dict d = result_dict(m, numerics::solve_linear_static(m, kind_of(solver)));
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
          "backend that actually ran the solve).");
  mod.def("solve_text", &solve_text, py::arg("text"), py::arg("solver") = "",
          py::arg("backend") = "",
          "Like solve() but takes the deck contents as a string.");

  mod.def("summary", &summary_file, py::arg("path"),
          "Parse an .inp deck without solving and return a dict describing it: "
          "num_nodes, num_elements, num_materials, materials (sorted names), and "
          "requested_solver ('direct'/'cg').");
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
