#include "calculixpp/numerics/buckling.hpp"

#include "calculixpp/core/results.hpp"
#include "calculixpp/fem/assembly.hpp"
#include "calculixpp/fem/stress.hpp"
#include "calculixpp/numerics/linear_static.hpp"

namespace cxpp::numerics {

BucklingReport solve_buckling(const Model& model, std::size_t num_modes) {
  // Step A — linear static prestress solve for the reference load, then recover the
  // per-element per-Gauss reference stress field feeding K_geo.
  const StaticFields prestress = solve_linear_static(model);
  const std::vector<std::vector<Voigt6>> gp_stress =
      fem::recover_gauss_stress(model, prestress.displacement);

  // Step B — assemble the unloaded stiffness K (SPD Cholesky anchor) and the geometric
  // stiffness K_geo(σ_ref) on the same free-DOF numbering, then extract the lowest
  // positive buckling factors of (K + λ K_geo) φ = 0.
  const fem::LinearSystem K = fem::assemble_linear_static(model);
  const fem::LinearSystem Kgeo = fem::assemble_geometric_stiffness(model, gp_stress);

  BucklingReport report;
  report.basis = extract_buckling_modes(K, Kgeo, num_modes);
  report.factors.reserve(report.basis.modes.size());
  for (const Mode& m : report.basis.modes) report.factors.push_back(m.eigenvalue);
  return report;
}

}  // namespace cxpp::numerics
