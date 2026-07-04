#include "calculixpp/numerics/high_cycle_fatigue.hpp"

#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>

#include "calculixpp/core/results.hpp"
#include "calculixpp/numerics/linear_static.hpp"

namespace cxpp::numerics {

namespace {

// von Mises equivalent of a Voigt6 stress (xx,yy,zz,xy,xz,yz, engineering shear).
Real von_mises(const Voigt6& s) {
  const Real a = s[0] - s[1], b = s[1] - s[2], c = s[2] - s[0];
  return std::sqrt(0.5 * (a * a + b * b + c * c) +
                   3.0 * (s[3] * s[3] + s[4] * s[4] + s[5] * s[5]));
}

// Scalar uniaxial-equivalent stress amplitude for the selected fatigue criterion. Signed
// von Mises multiplies the (non-negative) von Mises equivalent by the sign of the
// hydrostatic (mean normal) stress so tension/compression are distinguished; plain von
// Mises returns the equivalent unchanged.
Real amplitude(const Voigt6& s, FatigueCriterion crit) {
  const Real vm = von_mises(s);
  if (crit == FatigueCriterion::VonMises) return vm;
  const Real mean = (s[0] + s[1] + s[2]) / 3.0;
  const Real sign = mean < 0.0 ? -1.0 : 1.0;
  return sign * vm;
}

// First material carrying a *FATIGUE (Basquin S-N) curve, or nullopt when none does. The
// stress-life slice applies a single S-N curve across the model (per-element resolution is
// a follow-on; see proposal). A degenerate curve (SNCurve::empty) does not count.
std::optional<SNCurve> find_sn_curve(const Model& model) {
  for (const auto& [name, mat] : model.materials)
    if (mat.sn_curve && !mat.sn_curve->empty()) return *mat.sn_curve;
  return std::nullopt;
}

}  // namespace

HcfReport evaluate_hcf(const Model& model) {
  // A valid preceding results source needs (a) an S-N curve to invert and (b) an elastic
  // model to recover a stress field from. Diagnose either missing source before solving.
  const std::optional<SNCurve> sn = find_sn_curve(model);
  if (!sn)
    throw std::runtime_error(
        "*HCF: no S-N (*FATIGUE) curve found on any material — cannot evaluate fatigue "
        "life");

  bool has_elastic = false;
  for (const auto& [name, mat] : model.materials)
    if (mat.effective_elastic()) has_elastic = true;
  if (!has_elastic)
    throw std::runtime_error(
        "*HCF: no elastic material to recover a stress field from — missing stress "
        "source");

  // Recover the preceding stress field (linear-static path). `StaticFields::stress` is the
  // averaged nodal stress, aligned with the mesh node indices.
  const StaticFields fields = solve_linear_static(model);
  const std::size_t nn = model.mesh.num_nodes();

  HcfReport report;
  report.life.resize(nn, 0.0);
  report.amplitude.resize(nn, 0.0);

  // Basquin S_a = a * N^b  ->  N = (S_a / a)^(1/b). With b < 0, a larger amplitude yields a
  // smaller N (shorter life). A zero/near-zero amplitude gives an effectively infinite
  // life (represented as +inf); the worst case is the largest amplitude / smallest life.
  const Real inv_b = 1.0 / sn->b;
  Real worst_amp = -std::numeric_limits<Real>::infinity();
  std::size_t worst_i = 0;
  for (std::size_t i = 0; i < nn; ++i) {
    const Real s_a = amplitude(fields.stress[i], model.hcf_criterion);
    report.amplitude[i] = s_a;
    const Real abs_a = std::fabs(s_a);
    report.life[i] = abs_a > 0.0 ? std::pow(abs_a / sn->a, inv_b)
                                 : std::numeric_limits<Real>::infinity();
    if (abs_a > worst_amp) {
      worst_amp = abs_a;
      worst_i = i;
    }
  }

  report.worst_index = worst_i;
  report.worst_node_id = model.mesh.nodes()[worst_i].id;
  report.worst_location = model.mesh.nodes()[worst_i].x;
  report.worst_amplitude = report.amplitude[worst_i];
  report.worst_life = report.life[worst_i];
  return report;
}

}  // namespace cxpp::numerics
