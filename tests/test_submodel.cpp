// Submodeling validation (spec: submodeling). A coarse GLOBAL bar is solved; a finer
// SUBMODEL is cut from a sub-region and driven on its full cut boundary by the global
// displacement interpolated at each boundary node. The submodel must reproduce the
// global displacement field inside the region to a tight relative-L2 tolerance.
#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "calculixpp/core/model.hpp"
#include "calculixpp/core/submodel.hpp"
#include "calculixpp/numerics/linear_static.hpp"
#include "calculixpp/numerics/submodel.hpp"
#include "check.hpp"

using namespace cxpp;

namespace {

constexpr Real kE = 210000.0, kNu = 0.3;

// Build a structured hex (C3D8) grid mesh over [0,Lx]x[0,Ly]x[0,Lz] with nx*ny*nz
// cells. Nodes are numbered 1..(nx+1)(ny+1)(nz+1) in (i fastest, then j, then k) order.
struct Grid {
  int nx, ny, nz;
  Real Lx, Ly, Lz;
  Index nid(int i, int j, int k) const {
    return 1 + i + (nx + 1) * (j + (ny + 1) * k);
  }
  Vec3 pos(int i, int j, int k) const {
    return {Lx * i / nx, Ly * j / ny, Lz * k / nz};
  }
};

void build_grid_mesh(Model& m, const Grid& g) {
  for (int k = 0; k <= g.nz; ++k)
    for (int j = 0; j <= g.ny; ++j)
      for (int i = 0; i <= g.nx; ++i) m.mesh.add_node(g.nid(i, j, k), g.pos(i, j, k));
  Index eid = 1;
  std::vector<Index> all;
  for (int k = 0; k < g.nz; ++k)
    for (int j = 0; j < g.ny; ++j)
      for (int i = 0; i < g.nx; ++i) {
        // C3D8 node order: bottom face CCW then top face CCW (Abaqus convention).
        const std::vector<Index> conn = {
            g.nid(i, j, k),     g.nid(i + 1, j, k),     g.nid(i + 1, j + 1, k),
            g.nid(i, j + 1, k), g.nid(i, j, k + 1),     g.nid(i + 1, j, k + 1),
            g.nid(i + 1, j + 1, k + 1), g.nid(i, j + 1, k + 1)};
        m.mesh.add_element(eid, ElementType::C3D8, conn);
        all.push_back(eid);
        ++eid;
      }
  m.mesh.add_elset("EALL", all);
  m.materials["MAT"] = Material{"MAT", ElasticIso{kE, kNu}, std::nullopt};
  m.sections.push_back(SolidSection{"EALL", "MAT"});
}

// Capture a solved model's mesh + displacement as a GlobalSolution.
GlobalSolution to_global(const Model& m, const StaticFields& f) {
  GlobalSolution g;
  for (const Node& n : m.mesh.nodes()) {
    g.node_ids.push_back(n.id);
    g.coords.push_back(n.x);
  }
  g.displacement = f.displacement;
  for (const Element& e : m.mesh.elements()) {
    g.elem_type.push_back(e.type);
    std::vector<Index> conn;
    for (Index nid : e.nodes) conn.push_back(m.mesh.node_index(nid));
    g.elem_conn.push_back(std::move(conn));
  }
  return g;
}

// Relative L2 error of the submodel displacement vs. the "truth" evaluated at each
// submodel node's location.
Real rel_l2(const Model& sub, const StaticFields& fs,
            const std::vector<Vec3>& truth) {
  Real num = 0, den = 0;
  for (std::size_t i = 0; i < sub.mesh.num_nodes(); ++i)
    for (int c = 0; c < 3; ++c) {
      const Real d = fs.displacement[i][static_cast<std::size_t>(c)] -
                     truth[i][static_cast<std::size_t>(c)];
      num += d * d;
      den += truth[i][static_cast<std::size_t>(c)] * truth[i][static_cast<std::size_t>(c)];
    }
  return std::sqrt(num / (den > 0 ? den : 1.0));
}

// A submodel node is on the cut boundary when it lies on the submodel's outer box
// surface — every such node is driven from the global solution.
bool on_box_surface(const Vec3& x, Real x0, Real x1, Real y0, Real y1, Real z0,
                    Real z1) {
  const Real t = 1e-9;
  return std::fabs(x[0] - x0) < t || std::fabs(x[0] - x1) < t ||
         std::fabs(x[1] - y0) < t || std::fabs(x[1] - y1) < t ||
         std::fabs(x[2] - z0) < t || std::fabs(x[2] - z1) < t;
}

// GLOBAL: a bar along x, fixed at x=0, pulled at x=L. The linear-elastic response is a
// smooth FE field; a sub-region cut from it and driven on its boundary must reproduce
// it. (Uniaxial tension gives a linear u field that any conforming mesh reproduces
// exactly, so the oracle is tight and independent of mesh density.)
void run_validation() {
  const Grid gg{6, 2, 2, 6.0, 2.0, 2.0};
  Model global;
  build_grid_mesh(global, gg);
  // Fix the x=0 face (all DOFs on the built-in end), pull the x=L face in +x.
  for (int k = 0; k <= gg.nz; ++k)
    for (int j = 0; j <= gg.ny; ++j) {
      for (int c = 1; c <= 3; ++c) global.spcs.push_back(Spc{gg.nid(0, j, k), c, 0.0});
      global.cloads.push_back(Cload{gg.nid(gg.nx, j, k), 1, 500.0});
    }
  const StaticFields gf = numerics::solve_linear_static(global);
  const GlobalSolution gsol = to_global(global, gf);

  // Interpolation sanity: at a global node location the interpolated displacement equals
  // that node's solved displacement.
  {
    const Index probe = global.mesh.node_index(gg.nid(3, 1, 1));
    const Vec3 X = global.mesh.nodes()[static_cast<std::size_t>(probe)].x;
    const Vec3 u = numerics::interpolate_global_displacement(gsol, X);
    for (int c = 0; c < 3; ++c)
      CX_NEAR(u[static_cast<std::size_t>(c)],
              gf.displacement[static_cast<std::size_t>(probe)][static_cast<std::size_t>(c)],
              1e-9);
  }

  // SUBMODEL: a finer mesh over the interior sub-box [2,4]x[0.5,1.5]x[0.5,1.5].
  const Real x0 = 2.0, x1 = 4.0, y0 = 0.5, y1 = 1.5, z0 = 0.5, z1 = 1.5;
  Model sub;
  // Offset the sub-grid to the sub-box; build with a shifted Grid then translate.
  const Grid sg{4, 3, 3, x1 - x0, y1 - y0, z1 - z0};
  build_grid_mesh(sub, sg);
  // Translate the submodel nodes into the global sub-box position by rebuilding coords.
  // (add_node already placed them at origin-based coords; re-add with the offset.)
  Model sub2;
  for (int k = 0; k <= sg.nz; ++k)
    for (int j = 0; j <= sg.ny; ++j)
      for (int i = 0; i <= sg.nx; ++i) {
        const Vec3 p = sg.pos(i, j, k);
        sub2.mesh.add_node(sg.nid(i, j, k), {p[0] + x0, p[1] + y0, p[2] + z0});
      }
  Index eid = 1;
  std::vector<Index> all;
  for (int k = 0; k < sg.nz; ++k)
    for (int j = 0; j < sg.ny; ++j)
      for (int i = 0; i < sg.nx; ++i) {
        const std::vector<Index> conn = {
            sg.nid(i, j, k),     sg.nid(i + 1, j, k),     sg.nid(i + 1, j + 1, k),
            sg.nid(i, j + 1, k), sg.nid(i, j, k + 1),     sg.nid(i + 1, j, k + 1),
            sg.nid(i + 1, j + 1, k + 1), sg.nid(i, j + 1, k + 1)};
        sub2.mesh.add_element(eid++, ElementType::C3D8, conn);
        all.push_back(eid - 1);
      }
  sub2.mesh.add_elset("EALL", all);
  sub2.mesh.add_nset("BND", {});  // filled below
  sub2.materials["MAT"] = Material{"MAT", ElasticIso{kE, kNu}, std::nullopt};
  sub2.sections.push_back(SolidSection{"EALL", "MAT"});
  // Declare the submodel and drive every boundary node's 3 DOFs.
  sub2.submodels.push_back(SubmodelSpec{"BND", "EALL"});
  for (const Node& n : sub2.mesh.nodes())
    if (on_box_surface(n.x, x0, x1, y0, y1, z0, z1))
      for (int c = 1; c <= 3; ++c) {
        Spc spc{n.id, c, 0.0};
        spc.driven = true;
        sub2.spcs.push_back(spc);
      }

  const StaticFields sf = numerics::solve_submodel(sub2, gsol);

  // Truth: the global displacement interpolated at each submodel node's location.
  std::vector<Vec3> truth(sub2.mesh.num_nodes());
  for (std::size_t i = 0; i < sub2.mesh.num_nodes(); ++i)
    truth[i] = numerics::interpolate_global_displacement(
        gsol, sub2.mesh.nodes()[i].x);
  const Real err = rel_l2(sub2, sf, truth);
  std::fprintf(stderr, "submodel rel-L2 error = %.3e\n", err);
  CX_CHECK(err < 1e-3);
}

// A driven node outside the global element set must throw.
void run_outside() {
  const Grid gg{2, 1, 1, 1.0, 1.0, 1.0};
  Model global;
  build_grid_mesh(global, gg);
  for (int j = 0; j <= gg.ny; ++j)
    for (int k = 0; k <= gg.nz; ++k)
      for (int c = 1; c <= 3; ++c) global.spcs.push_back(Spc{gg.nid(0, j, k), c, 0.0});
  global.cloads.push_back(Cload{gg.nid(gg.nx, 0, 0), 1, 100.0});
  const StaticFields gf = numerics::solve_linear_static(global);
  const GlobalSolution gsol = to_global(global, gf);
  bool threw = false;
  try {
    numerics::interpolate_global_displacement(gsol, {10.0, 10.0, 10.0});
  } catch (const std::exception&) {
    threw = true;
  }
  CX_CHECK(threw);
}

}  // namespace

int main() {
  run_validation();
  run_outside();
  CX_MAIN_RETURN();
}
