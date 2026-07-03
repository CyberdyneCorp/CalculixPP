// ccxpp — command-line runner for the CalculiX++ linear-static solver.
//   ccxpp <model.inp> [-o <basename>] [--solver direct|cg]
// Parses an Abaqus-style deck, solves K u = f, and writes <basename>.frd/.dat.
#include <cmath>
#include <cstdio>
#include <optional>
#include <string>

#include "calculixpp/io/inp_parser.hpp"
#include "calculixpp/io/results_writer.hpp"
#include "calculixpp/numerics/heat_transfer.hpp"
#include "calculixpp/numerics/linear_static.hpp"
#include "calculixpp/numerics/multistep.hpp"
#include "calculixpp/numerics/nonlinear_static.hpp"

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
  static const char* kUsage =
      "usage: ccxpp <model.inp> [-o <basename>] [--solver direct|cg] "
      "[--nonlinear] [--line-search]\n";
  std::string input, out;
  std::optional<numerics::SolverKind> forced;  // empty -> deck SOLVER= / Auto
  bool nonlinear = false, line_search = false;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "-o" && i + 1 < argc) {
      out = argv[++i];
    } else if (a == "--solver" && i + 1 < argc) {
      const std::string s = argv[++i];
      forced = (s == "cg") ? numerics::SolverKind::CG : numerics::SolverKind::Direct;
    } else if (a == "--nonlinear") {
      nonlinear = true;
    } else if (a == "--line-search") {
      line_search = true;
    } else if (!a.empty() && a[0] != '-') {
      input = a;
    } else {
      std::fprintf(stderr, "%s", kUsage);
      return 2;
    }
  }
  if (input.empty()) {
    std::fprintf(stderr, "%s", kUsage);
    return 2;
  }
  if (out.empty()) out = stem(input);

  try {
    // Multi-step mechanical decks (>1 *STEP) run the step-loop driver and write the
    // final step's fields. A single-*STEP deck falls through to the byte-identical
    // single-step path below (parse_inp_steps returns one model there). Thermal/coupled
    // multi-step is out of this slice; those decks are handled single-step as before.
    {
      const std::vector<Model> steps = cxpp::io::parse_inp_steps_file(input);
      if (steps.size() > 1 && steps.front().procedure == Procedure::Static &&
          !nonlinear && !steps.front().has_nonlinear_material()) {
        const Model& last = steps.back();
        const StaticFields f = numerics::solve_multistep_static(steps, forced);
        Real umax = 0.0, svm_max = 0.0;
        for (std::size_t i = 0; i < last.mesh.num_nodes(); ++i) {
          const Vec3& u = f.displacement[i];
          umax = std::max(umax, std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
          svm_max = std::max(svm_max, von_mises(f.stress[i]));
        }
        cxpp::io::write_frd(out + ".frd", last, f);
        cxpp::io::write_dat(out + ".dat", last, f);
        std::printf("CalculiX++  %s  (multi-step: %zu steps)\n", input.c_str(),
                    steps.size());
        std::printf("  nodes=%zu  elements=%zu\n", last.mesh.num_nodes(),
                    last.mesh.num_elements());
        std::printf("  max |u|        = %.6g  (final step)\n", umax);
        std::printf("  max von Mises  = %.6g  (final step)\n", svm_max);
        std::printf("  wrote %s.frd, %s.dat\n", out.c_str(), out.c_str());
        return 0;
      }
    }

    const Model model = cxpp::io::parse_inp_file(input);

    // A *COUPLED TEMPERATURE-DISPLACEMENT deck solves the thermal field, applies the
    // resulting thermal strain, then solves the mechanical field (one-way coupling).
    // The written .frd/.dat carry the mechanical (displacement/stress) result.
    if (model.procedure == Procedure::Coupled) {
      const CoupledFields cf = numerics::solve_coupled(model, forced);
      Real umax = 0.0, svm_max = 0.0;
      for (std::size_t i = 0; i < model.mesh.num_nodes(); ++i) {
        const Vec3& u = cf.mechanical.displacement[i];
        umax = std::max(umax, std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
        svm_max = std::max(svm_max, von_mises(cf.mechanical.stress[i]));
      }
      cxpp::io::write_frd(out + ".frd", model, cf.mechanical);
      cxpp::io::write_dat(out + ".dat", model, cf.mechanical);
      std::printf("CalculiX++  %s  (coupled temperature-displacement)\n",
                  input.c_str());
      std::printf("  nodes=%zu  elements=%zu\n", model.mesh.num_nodes(),
                  model.mesh.num_elements());
      std::printf("  max |u|        = %.6g\n", umax);
      std::printf("  max von Mises  = %.6g\n", svm_max);
      std::printf("  wrote %s.frd, %s.dat\n", out.c_str(), out.c_str());
      return 0;
    }

    // Auto-dispatch on the parsed procedure: a *HEAT TRANSFER deck (steady or
    // transient) solves the scalar temperature field; everything else takes the
    // mechanical path.
    if (model.procedure == Procedure::HeatTransferSteady ||
        model.procedure == Procedure::HeatTransferTransient) {
      const bool transient = model.procedure == Procedure::HeatTransferTransient;
      const ThermalFields t = numerics::solve_heat_transfer(model, forced);
      Real tmin = 0.0, tmax = 0.0;
      if (!t.temperature.empty()) tmin = tmax = t.temperature[0];
      for (const Real v : t.temperature) {
        tmin = std::min(tmin, v);
        tmax = std::max(tmax, v);
      }
      cxpp::io::write_frd(out + ".frd", model, t);
      cxpp::io::write_dat(out + ".dat", model, t);
      std::printf("CalculiX++  %s  (heat transfer, %s)\n", input.c_str(),
                  transient ? "transient" : "steady state");
      std::printf("  nodes=%zu  elements=%zu\n", model.mesh.num_nodes(),
                  model.mesh.num_elements());
      std::printf("  temperature range = [%.6g, %.6g]\n", tmin, tmax);
      std::printf("  wrote %s.frd, %s.dat\n", out.c_str(), out.c_str());
      return 0;
    }

    // --solver overrides; otherwise honor the deck's SOLVER= (Auto by default).
    // --nonlinear routes through the Newton-Raphson driver (identical results on a
    // linear model); the default path stays the linear solve.
    StaticFields f;
    // A *CONTACT PAIR deck is a nonlinear constraint problem: route it to the
    // Newton driver, which assembles the penalty contact contribution into the
    // tangent + residual. --nonlinear or a nonlinear material also route here.
    if (nonlinear || model.has_nonlinear_material() || model.has_contact()) {
      numerics::NonlinearOptions opts;
      opts.line_search = line_search;
      opts.forced = forced;
      f = numerics::solve_nonlinear_static(model, opts);
    } else {
      f = numerics::solve_linear_static(model, forced);
    }

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
