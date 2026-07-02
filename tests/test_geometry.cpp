// Geometry stub tests (spec: mesh-processing). CyberCadKernel is unavailable in
// Phase 1, so every CadMesher operation must throw the documented message. This
// documents the bypass path: the real Phase-1 mesh comes from deck import.
#include <stdexcept>
#include <string>

#include "calculixpp/core/mesh.hpp"
#include "calculixpp/geometry/cad_mesher.hpp"
#include "check.hpp"

using namespace cxpp;
using namespace cxpp::geometry;

namespace {

// Run `fn`; verify it throws std::runtime_error carrying kCadUnavailableMessage.
template <typename Fn>
void expect_unavailable(Fn&& fn) {
  bool threw = false;
  try {
    fn();
  } catch (const std::runtime_error& e) {
    threw = true;
    CX_CHECK(std::string(e.what()) == kCadUnavailableMessage);
  }
  CX_CHECK(threw);
}

void test_all_ops_throw_documented_message() {
  const CadMesher m;
  const BoundaryRep brep{};
  const SurfaceMesh surf{};
  Mesh mesh;

  expect_unavailable([&] { (void)m.import("part.step", CadFormat::Step); });
  expect_unavailable([&] { (void)m.heal(brep, HealingOptions{}); });
  expect_unavailable([&] { (void)m.triangulate(brep, TriangulationOptions{}); });
  expect_unavailable([&] { (void)m.tet_mesh(surf, VolumeMeshOptions{}); });
  expect_unavailable([&] { (void)m.quality(mesh); });
  expect_unavailable([&] { m.map_to_model(mesh, mesh); });
}

}  // namespace

int main() {
  test_all_ops_throw_documented_message();
  if (cxtest::g_failures == 0) std::printf("test_geometry: OK\n");
  CX_MAIN_RETURN();
}
