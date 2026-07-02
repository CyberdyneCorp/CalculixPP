// ccxpp — command-line runner for the CalculiX++ linear-static solver.
//   ccxpp <model.inp> [-o <basename>] [--solver direct|cg]
// Parses an Abaqus-style deck, solves K u = f, and writes <basename>.frd/.dat.
#include <cmath>
#include <cstdio>
#include <optional>
#include <string>

#include "calculixpp/io/inp_parser.hpp"
#include "calculixpp/io/results_writer.hpp"
#include "calculixpp/numerics/linear_static.hpp"

namespace {

using namespace cxpp;

Real von_mises(const Voigt6& s) {
  const Real a = s[0] - s[1], b = s[1] - s[2], c = s[2] - s[0];
  return std::sqrt(0.5 * (a * a + b * b + c * c) +
                   3.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]));
}

std::string stem(const std::string& path) {
  const auto slash = path.find_last_of("/\\");
  const std::string base = slash == std::string::npos ? path : path.substr(slash + 1);
  const auto dot = base.find_last_of('.');
  return dot == std::string::npos ? base : base.substr(0, dot);
}

}  // namespace

int main(int argc, char** argv) {
  std::string input, out;
  std::optional<numerics::SolverKind> forced;  // empty -> deck SOLVER= / Auto
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "-o" && i + 1 < argc) {
      out = argv[++i];
    } else if (a == "--solver" && i + 1 < argc) {
      const std::string s = argv[++i];
      forced = (s == "cg") ? numerics::SolverKind::CG : numerics::SolverKind::Direct;
    } else if (!a.empty() && a[0] != '-') {
      input = a;
    } else {
      std::fprintf(stderr, "usage: ccxpp <model.inp> [-o <basename>] [--solver direct|cg]\n");
      return 2;
    }
  }
  if (input.empty()) {
    std::fprintf(stderr, "usage: ccxpp <model.inp> [-o <basename>] [--solver direct|cg]\n");
    return 2;
  }
  if (out.empty()) out = stem(input);

  try {
    const Model model = cxpp::io::parse_inp_file(input);
    // --solver overrides; otherwise honor the deck's SOLVER= (Auto by default).
    const StaticFields f = numerics::solve_linear_static(model, forced);

    Real umax = 0.0, svm_max = 0.0;
    for (std::size_t i = 0; i < model.mesh.num_nodes(); ++i) {
      const Vec3& u = f.displacement[i];
      umax = std::max(umax, std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
      svm_max = std::max(svm_max, von_mises(f.stress[i]));
    }

    cxpp::io::write_frd(out + ".frd", model, f);
    cxpp::io::write_dat(out + ".dat", model, f);

    std::printf("CalculiX++  %s\n", input.c_str());
    std::printf("  nodes=%zu  elements=%zu\n", model.mesh.num_nodes(),
                model.mesh.num_elements());
    std::printf("  max |u|        = %.6g\n", umax);
    std::printf("  max von Mises  = %.6g\n", svm_max);
    std::printf("  wrote %s.frd, %s.dat\n", out.c_str(), out.c_str());
    return 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "ccxpp: error: %s\n", e.what());
    return 1;
  }
}
