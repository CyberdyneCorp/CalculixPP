#pragma once
#include <string>
#include <vector>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/types.hpp"
#include "calculixpp/fem/assembly.hpp"

// Scalar thermal (heat-conduction) assembly — a PARALLEL path to the mechanical
// assembly with ONE DOF per node (the temperature). It reuses the same element
// shape functions / Gauss rules (via element_conduction), the same COO->SciPP solve
// (numerics::solve_reduced), and the same LinearSystem container, but numbers a
// single temperature DOF per node and eliminates prescribed temperatures directly
// (Dirichlet static condensation, mirroring the mechanical SPC path but scalar).
// The mechanical path is untouched. (spec: heat-transfer-analysis — steady conduction.)
namespace cxpp::fem {

// Assemble the reduced steady conduction system Kt_ff T_f = q_f - Kt_fc T_c for a
// heat-transfer step: the free/free conductivity in COO form and the reduced thermal
// rhs (concentrated + distributed fluxes minus the eliminated prescribed-temperature
// columns). The returned LinearSystem uses a SCALAR DOF map: dof_eq[node_index] is
// the free equation (>=0) or -1 if the node's temperature is prescribed; prescribed[
// node_index] holds the fixed temperature. transform/prescribed_amp are unused
// (no MPCs on the thermal path). Film (*FILM) contributions are included: they add
// h ∫N_a N_b dA to the matrix and h T_sink ∫N_a dA to the rhs. Radiation is NOT
// included here (it is nonlinear; the transient/steady drivers add it iteratively).
// (spec: heat-transfer-analysis.)
LinearSystem assemble_conduction(const Model& model);

// Scalar thermal DOF map (temperature per node): fills sys.dof_eq (free equation or
// -1 if prescribed), sys.prescribed (fixed temperature per node), sys.n_free, and a
// zeroed sys.rhs. Shared by every thermal driver so they number DOFs identically.
void build_thermal_dof_map(const Model& model, LinearSystem& sys);

// The thermal operators on the FULL node numbering (before Dirichlet reduction),
// used by the transient driver. `conductivity` and `capacitance` are sparse
// triplet lists (COO) over node indices; `flux` is the applied source per node
// (concentrated + distributed flux + film source h T_sink ∫N). Film matrix terms
// h ∫N_a N_b are folded into `conductivity`. Radiation is added by the driver
// per iteration (nonlinear), not here. Assembled at step fraction `lambda`.
struct FullThermalSystem {
  std::vector<Index> k_rows, k_cols;  // conductivity + film COO
  std::vector<Real> k_vals;
  std::vector<Index> c_rows, c_cols;  // capacitance COO
  std::vector<Real> c_vals;
  std::vector<Real> flux;             // source per node, size num_nodes
};
FullThermalSystem assemble_full_thermal(const Model& model, Real lambda);

}  // namespace cxpp::fem
