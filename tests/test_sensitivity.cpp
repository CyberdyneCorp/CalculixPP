// Design-sensitivity core (add-design-optimization). The authoritative correctness
// gate is a CENTRAL FINITE-DIFFERENCE gradient check: the adjoint dObjective/dx must
// match (O(+h) − O(−h)) / 2h from a re-solve at perturbed coordinates, to < 1e-4
// relative, on a small cantilever beam — for BOTH a compliance response (self-adjoint)
// and a nodal-displacement response (general linear adjoint). Also checks the parser
// wiring (design variables / response) and that deferred cards are rejected.
#include <cmath>
#include <string>
#include <vector>

#include "calculixpp/fem/assembly.hpp"
#include "calculixpp/io/inp_parser.hpp"
#include "calculixpp/numerics/linear_static.hpp"
#include "calculixpp/numerics/sensitivity.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

// A short cantilever built from two C3D8 hex bricks along x, fixed at x=0, tip-loaded
// in −y. The two mid/tip faces (x=1,2) give interior + tip design nodes to perturb.
// Node numbering: face k (k=0,1,2) at x=k has nodes 4k+1..4k+4 (y,z corners).
std::string beam_deck(const std::string& response_card) {
  return std::string(R"(
*NODE, NSET=NALL
 1, 0.0, 0.0, 0.0
 2, 0.0, 1.0, 0.0
 3, 0.0, 1.0, 1.0
 4, 0.0, 0.0, 1.0
 5, 1.0, 0.0, 0.0
 6, 1.0, 1.0, 0.0
 7, 1.0, 1.0, 1.0
 8, 1.0, 0.0, 1.0
 9, 2.0, 0.0, 0.0
 10, 2.0, 1.0, 0.0
 11, 2.0, 1.0, 1.0
 12, 2.0, 0.0, 1.0
*ELEMENT, TYPE=C3D8, ELSET=EALL
 1, 1, 5, 6, 2, 4, 8, 7, 3
 2, 5, 9, 10, 6, 8, 12, 11, 7
*NSET, NSET=FIX
 1, 2, 3, 4
*NSET, NSET=DVNODES
 6, 7, 10, 11
*MATERIAL, NAME=STEEL
*ELASTIC
 210000.0, 0.3
*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL
*BOUNDARY
 FIX, 1, 3, 0.0
*STEP
*STATIC
*CLOAD
 10, 2, -100.0
 11, 2, -100.0
*END STEP
*STEP
*SENSITIVITY
*DESIGN VARIABLES, TYPE=COORDINATE
 DVNODES
)") + response_card + "*END STEP\n";
}

// Evaluate the objective the SAME way solve_sensitivity defines it, at a model whose
// design coordinate (var) is shifted by `delta`. Compliance = b·u; Displacement = u[eq].
Real objective_at(const Model& base, const DesignVariable& var, Real delta,
                  const DesignResponse& resp) {
  Model m = base;
  const Index ni = m.mesh.node_index(var.node_id);
  if (delta != 0.0) m.mesh.perturb_node_coord(ni, var.comp - 1, delta);
  const fem::LinearSystem sys = fem::assemble_linear_static(m);
  const std::vector<Real> u = numerics::solve_reduced(sys, numerics::SolverKind::Direct);
  if (resp.kind == DesignResponse::Kind::Compliance) {
    Real g = 0.0;
    for (std::size_t i = 0; i < u.size(); ++i) g += sys.rhs[i] * u[i];
    return g;
  }
  const Index rni = m.mesh.node_index(resp.node_id);
  const Index eq = sys.dof_eq[static_cast<std::size_t>(rni) * kDofsPerNode +
                              static_cast<std::size_t>(resp.comp - 1)];
  return u[static_cast<std::size_t>(eq)];
}

// Central FD of one response over one design variable at step `hfd`.
Real fd_gradient(const Model& m, const DesignVariable& var, const DesignResponse& resp,
                 Real hfd) {
  return (objective_at(m, var, hfd, resp) - objective_at(m, var, -hfd, resp)) /
         (2.0 * hfd);
}

// The authoritative gate: adjoint dObjective/dx vs central FD, < 1e-4 relative.
void check_fd_gradient(const std::string& response_card, const char* label) {
  const Model m = io::parse_inp(beam_deck(response_card));
  CX_CHECK(m.procedure == Procedure::Sensitivity);
  CX_CHECK(!m.design_variables.empty());
  CX_CHECK(!m.design_responses.empty());

  const numerics::SensitivityReport rep = numerics::solve_sensitivity(m);
  CX_CHECK(rep.responses.size() == 1);
  const numerics::SensitivityResult& r = rep.responses.front();
  CX_CHECK(r.dgdx.size() == m.design_variables.size());

  const Real hfd = 1e-5;  // FD step (independent of the driver's semi-analytic h)
  Real max_rel = 0.0;
  bool any_significant = false;
  for (std::size_t v = 0; v < m.design_variables.size(); ++v) {
    const Real fd = fd_gradient(m, m.design_variables[v], m.design_responses.front(),
                                hfd);
    const Real adj = r.dgdx[v];
    const Real scale = std::max(std::fabs(fd), std::fabs(adj));
    if (scale < 1e-8) continue;  // both ~0: no meaningful relative error
    any_significant = true;
    const Real rel = std::fabs(fd - adj) / scale;
    max_rel = std::max(max_rel, rel);
    if (rel > 1e-4)
      std::fprintf(stderr,
                   "  [%s] var %zu (node %lld comp %d): adjoint=%.8e fd=%.8e rel=%.2e\n",
                   label, v, static_cast<long long>(m.design_variables[v].node_id),
                   m.design_variables[v].comp, adj, fd, rel);
  }
  std::fprintf(stderr, "[%s] max relative gradient error = %.3e\n", label, max_rel);
  CX_CHECK(any_significant);
  CX_CHECK(max_rel < 1e-4);
}

// Deferred cards (*FILTER, *FEASIBLE DIRECTION, TYPE=ORIENTATION) must be rejected,
// never silently ignored / faked.
void test_deferred_rejected() {
  bool threw = false;
  try {
    io::parse_inp(
        "*NODE, NSET=N\n 1,0,0,0\n*STEP\n*SENSITIVITY\n*FILTER\n 1.0\n*END STEP\n");
  } catch (const std::exception&) {
    threw = true;
  }
  CX_CHECK(threw);

  threw = false;
  try {
    io::parse_inp(
        "*NODE, NSET=N\n 1,0,0,0\n*STEP\n*SENSITIVITY\n"
        "*DESIGN VARIABLES, TYPE=ORIENTATION\n N\n*END STEP\n");
  } catch (const std::exception&) {
    threw = true;
  }
  CX_CHECK(threw);
}

}  // namespace

int main() {
  // Compliance response (self-adjoint): STRAIN ENERGY over the whole model.
  check_fd_gradient("*DESIGN RESPONSE, NAME=COMPL\nSTRAIN ENERGY, EALL\n", "compliance");
  // Nodal-displacement response (general linear adjoint): tip node 10, y-DOF.
  check_fd_gradient("*DESIGN RESPONSE, NAME=TIPDISP\nU, 10, 2\n", "displacement");
  test_deferred_rejected();
  CX_MAIN_RETURN();
}
