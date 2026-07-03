// Newton-Raphson driver acceptance (spec: nonlinear-solution-control 1.5 GATE).
// For a linear model, solve_nonlinear_static must reproduce solve_linear_static to
// rel-L2 < 1e-10, converging in one increment / ~2 iterations. Also exercises the
// incrementation engine (multi-increment ramp) and the line-search flag: both must
// land on the same linear solution.
#include <cmath>
#include <vector>

#include "calculixpp/numerics/linear_static.hpp"
#include "calculixpp/numerics/nonlinear_static.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

Model make_single_tet() {
  Model m;
  m.mesh.add_node(1, {0, 0, 0});
  m.mesh.add_node(2, {1, 0, 0});
  m.mesh.add_node(3, {0, 1, 0});
  m.mesh.add_node(4, {0, 0, 1});
  m.mesh.add_element(100, ElementType::C3D4, {1, 2, 3, 4});
  m.mesh.add_elset("EALL", {100});
  m.materials["EL"] = Material{"EL", ElasticIso{210000.0, 0.3}, std::nullopt};
  m.sections.push_back(SolidSection{"EALL", "EL"});
  for (Index nd : {1, 2, 3})
    for (int c = 1; c <= 3; ++c) m.spcs.push_back(Spc{nd, c, 0.0});
  m.cloads.push_back(Cload{4, 3, -1000.0});
  return m;
}

Real rel_l2(const std::vector<Vec3>& a, const std::vector<Vec3>& b) {
  Real num = 0.0, den = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i)
    for (int c = 0; c < 3; ++c) {
      const Real d = a[i][static_cast<std::size_t>(c)] - b[i][static_cast<std::size_t>(c)];
      num += d * d;
      den += b[i][static_cast<std::size_t>(c)] * b[i][static_cast<std::size_t>(c)];
    }
  return den > 0 ? std::sqrt(num / den) : std::sqrt(num);
}

void check_reproduces_linear(const Model& m, const numerics::NonlinearOptions& opts) {
  const StaticFields lin = numerics::solve_linear_static(m);
  numerics::NonlinearReport rep;
  const StaticFields nl = numerics::solve_nonlinear_static(m, opts, &rep);
  CX_CHECK(nl.displacement.size() == lin.displacement.size());
  CX_CHECK(rep.converged);
  CX_CHECK(rep.final_load_factor > 1.0 - 1e-12);
  CX_CHECK(rel_l2(nl.displacement, lin.displacement) < 1e-10);
}

}  // namespace

int main() {
  const Model m = make_single_tet();

  // GATE: default single increment reproduces the linear solve exactly, in one
  // increment and 2 iterations (u=0 -> r=f_ext -> du=u_lin; then r=0).
  {
    numerics::NonlinearReport rep;
    const StaticFields lin = numerics::solve_linear_static(m);
    const StaticFields nl = numerics::solve_nonlinear_static(m, {}, &rep);
    CX_CHECK(rel_l2(nl.displacement, lin.displacement) < 1e-10);
    CX_CHECK(rep.increments == 1);
    CX_CHECK(rep.iterations == 2);
    CX_CHECK(rep.cutbacks == 0);
    CX_CHECK(rep.converged);
  }

  // Line search on: same linear solution.
  {
    numerics::NonlinearOptions opts;
    opts.line_search = true;
    check_reproduces_linear(m, opts);
  }

  // Automatic incrementation over several increments (initial 0.25) still lands on
  // the exact linear solution at load factor 1.
  {
    Model mm = m;
    mm.increment.initial = 0.25;
    mm.increment.max = 0.25;
    mm.increment.grow = 1.0;  // keep the step fixed to force multiple increments
    numerics::NonlinearReport rep;
    const StaticFields lin = numerics::solve_linear_static(mm);
    const StaticFields nl = numerics::solve_nonlinear_static(mm, {}, &rep);
    CX_CHECK(rep.increments >= 4);
    CX_CHECK(rep.converged);
    CX_CHECK(rel_l2(nl.displacement, lin.displacement) < 1e-10);
  }

  // DIRECT fixed increment: exactly 1/initial increments, no resizing.
  {
    Model mm = m;
    mm.increment.direct = true;
    mm.increment.initial = 0.5;
    numerics::NonlinearReport rep;
    const StaticFields lin = numerics::solve_linear_static(mm);
    const StaticFields nl = numerics::solve_nonlinear_static(mm, {}, &rep);
    CX_CHECK(rep.increments == 2);
    CX_CHECK(rep.cutbacks == 0);
    CX_CHECK(rel_l2(nl.displacement, lin.displacement) < 1e-10);
  }

  // Prescribed non-zero BC (enforced displacement) also reproduces the linear
  // solve: node 4 pulled to w=0.01, load removed.
  {
    Model mm = make_single_tet();
    mm.cloads.clear();
    mm.spcs.push_back(Spc{4, 3, 0.01});
    check_reproduces_linear(mm, {});
  }

  if (cxtest::g_failures == 0) std::printf("test_nonlinear: OK\n");
  CX_MAIN_RETURN();
}
