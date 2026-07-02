#pragma once
#include <string>
#include <vector>

#include "calculixpp/core/mesh.hpp"
#include "calculixpp/core/types.hpp"

// Geometry / mesh-processing module (spec: mesh-processing).
//
// PHASE 1 STATUS: interface stub only. The full pipeline — CAD import
// (STEP/IGES/BREP) -> boundary representation -> healing -> surface
// triangulation -> tetrahedral volume meshing -> mapping into the FE `Mesh` —
// is provided by CyberCadKernel, which is NOT available in the Phase-1 build.
// Every operation below is declared to pin the intended interface but is
// implemented as a stub that throws std::runtime_error with the message
// `kCadUnavailableMessage`. The working Phase-1 path is deck import: load the
// mesh directly from a `*NODE`/`*ELEMENT` deck via the input parser (see the
// mesh-processing "Import of existing analysis meshes" requirement), which
// bypasses CyberCadKernel entirely.
//
// This header carries NO CyberCadKernel dependency; it compiles and links into
// the dependency-free core so the intended API is visible to callers, tests,
// and (later) the Python bindings before the kernel is integrated.
namespace cxpp::geometry {

// Documented message thrown by every stub in this module. Callers and tests
// key off this exact text; keep it stable.
inline constexpr const char* kCadUnavailableMessage =
    "CyberCadKernel geometry meshing is not available in Phase 1; import mesh "
    "from the deck instead";

// Supported CAD exchange formats for import.
enum class CadFormat { Step, Iges, Brep };

// Requested tetrahedral element order for volume meshing.
enum class MeshOrder { Linear, Quadratic };  // C3D4, C3D10

// User-controllable healing tolerances (spec: Geometry healing and validation).
struct HealingOptions {
  Real gap_tolerance{1e-6};    // close gaps between faces below this size
  Real min_edge_length{1e-6};  // collapse edges shorter than this
  Real sliver_tolerance{1e-6}; // remove sliver faces below this
};

// User-controllable surface triangulation controls (spec: Surface triangulation).
struct TriangulationOptions {
  Real max_edge_length{0.0};    // 0 => kernel default; max triangle edge length
  Real chordal_deviation{0.0};  // 0 => kernel default; max chordal tolerance
  Real angular_deviation{0.0};  // 0 => kernel default; max angular tolerance (rad)
};

// User-controllable tetrahedral meshing controls (spec: Tetrahedral volume meshing).
struct VolumeMeshOptions {
  MeshOrder order{MeshOrder::Quadratic};
  Real target_element_size{0.0};  // 0 => kernel default
  Real grading{1.0};              // size grading factor
  Real min_scaled_jacobian{0.0};  // reject/refine below this quality
};

// Per-mesh quality summary (spec: Mesh quality metrics and reporting).
struct QualityReport {
  Real min_dihedral_angle{0.0};
  Real max_dihedral_angle{0.0};
  Real min_scaled_jacobian{0.0};
  Real mean_scaled_jacobian{0.0};
  Real max_aspect_ratio{0.0};
  std::size_t elements_below_threshold{0};
  std::vector<Index> flagged_elements;  // ids below the quality threshold
};

// Opaque boundary representation (B-rep) handle produced by import. The concrete
// shell/face/edge/vertex data lives inside CyberCadKernel; Phase 1 exposes only
// the reported summary so the interface is complete without the kernel type.
struct BoundaryRep {
  CadFormat format{CadFormat::Step};
  std::string units;              // reported geometric units
  std::size_t face_count{0};
  std::size_t edge_count{0};
  std::size_t vertex_count{0};
  std::size_t shell_count{0};
  std::vector<std::string> unsupported_entities;  // skipped, reported entities
};

// A triangulated surface mesh handle (spec: Surface triangulation). Watertight
// when the source B-rep is a valid closed shell.
struct SurfaceMesh {
  std::size_t triangle_count{0};
  std::size_t node_count{0};
  bool watertight{false};
  std::vector<std::string> free_edges;  // identified when open/non-manifold
};

// CyberCadKernel-backed geometry mesher. All members are stubs in Phase 1 and
// throw std::runtime_error(kCadUnavailableMessage). The real implementations
// land in a later phase behind the CyberCadKernel link.
class CadMesher {
 public:
  CadMesher() = default;

  // Import CAD geometry into an in-memory B-rep (spec: CAD geometry import).
  [[noreturn]] BoundaryRep import(const std::string& path, CadFormat format) const;

  // Validate and repair a B-rep before meshing (spec: Geometry healing).
  [[noreturn]] BoundaryRep heal(const BoundaryRep& brep,
                                const HealingOptions& options) const;

  // Tessellate B-rep faces into a triangle surface mesh (spec: Surface triangulation).
  [[noreturn]] SurfaceMesh triangulate(const BoundaryRep& brep,
                                       const TriangulationOptions& options) const;

  // Generate a tetrahedral volume mesh from a closed surface
  // (spec: Tetrahedral volume meshing).
  [[noreturn]] Mesh tet_mesh(const SurfaceMesh& surface,
                             const VolumeMeshOptions& options) const;

  // Compute per-element quality metrics for a generated mesh
  // (spec: Mesh quality metrics and reporting).
  [[noreturn]] QualityReport quality(const Mesh& mesh) const;

  // Map a generated mesh into the solver FE model (spec: Mapping into the FE model).
  [[noreturn]] void map_to_model(const Mesh& source, Mesh& out) const;
};

}  // namespace cxpp::geometry
