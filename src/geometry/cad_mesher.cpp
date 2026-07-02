// CyberCadKernel geometry / mesh-processing stubs (spec: mesh-processing).
//
// PHASE 1: CyberCadKernel is not available, so every entry point throws with
// the documented message. The working Phase-1 mesh path is deck import via the
// input parser (*NODE/*ELEMENT), which needs none of this. The real bodies land
// in a later phase once the kernel is linked.
#include "calculixpp/geometry/cad_mesher.hpp"

#include <stdexcept>

namespace cxpp::geometry {

namespace {
[[noreturn]] void unavailable() {
  throw std::runtime_error(kCadUnavailableMessage);
}
}  // namespace

BoundaryRep CadMesher::import(const std::string&, CadFormat) const { unavailable(); }

BoundaryRep CadMesher::heal(const BoundaryRep&, const HealingOptions&) const {
  unavailable();
}

SurfaceMesh CadMesher::triangulate(const BoundaryRep&,
                                   const TriangulationOptions&) const {
  unavailable();
}

Mesh CadMesher::tet_mesh(const SurfaceMesh&, const VolumeMeshOptions&) const {
  unavailable();
}

QualityReport CadMesher::quality(const Mesh&) const { unavailable(); }

void CadMesher::map_to_model(const Mesh&, Mesh&) const { unavailable(); }

}  // namespace cxpp::geometry
