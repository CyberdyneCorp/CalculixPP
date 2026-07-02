// Python bindings for CalculiX++ (spec: python-bindings).
// Exposes deck parsing and the linear-static solve, returning results as NumPy
// arrays for scripting and regression testing. C++ exceptions (incl. ParseError)
// propagate to Python as RuntimeError via pybind11's std::exception translation.
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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

numerics::SolverKind kind_of(const std::string& s) {
  return s == "cg" ? numerics::SolverKind::CG : numerics::SolverKind::Direct;
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
  d["reaction"] = vec3_to_array(f.reaction);
  d["num_nodes"] = n;
  d["num_elements"] = m.mesh.num_elements();
  return d;
}

py::dict solve_file(const std::string& path, const std::string& solver) {
  const Model m = io::parse_inp_file(path);
  return result_dict(m, numerics::solve_linear_static(m, kind_of(solver)));
}

py::dict solve_text(const std::string& text, const std::string& solver) {
  const Model m = io::parse_inp(text);
  return result_dict(m, numerics::solve_linear_static(m, kind_of(solver)));
}

}  // namespace

PYBIND11_MODULE(calculixpp, mod) {
  mod.doc() = "CalculiX++ — modern C++20 finite element solver (Python bindings)";
  mod.def("solve", &solve_file, py::arg("path"), py::arg("solver") = "direct",
          "Parse an .inp deck and solve the linear-static step. Returns a dict of "
          "NumPy arrays: node_ids, node_coords, displacement (Nx3), stress (Nx6, "
          "xx/yy/zz/xy/xz/yz), reaction (Nx3), num_nodes, num_elements.");
  mod.def("solve_text", &solve_text, py::arg("text"), py::arg("solver") = "direct",
          "Like solve() but takes the deck contents as a string.");
}
