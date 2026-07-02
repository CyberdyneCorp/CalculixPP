#pragma once
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "calculixpp/core/types.hpp"

// ComputeBackend abstraction (spec: compute-backend, linear-algebra-and-solvers).
//
// A pluggable surface over the numeric hot path. Phase 1 exposes only the sparse
// SPD solve used by the linear-static slice; assembly scatter / SpMV land as the
// backend matures. The CPU reference backend (src/compute/cpu_backend.cpp) wraps
// the SciPP sparse path (scipp::sparse from_coo -> spsolve/cg) and is the default,
// always-available implementation. GPU backends (CUDA/OpenCL/Metal) are future
// work and never required: absence of a GPU toolkit degrades to CPU, never to a
// build or run error.
namespace cxpp::compute {

// Which sparse solver the backend should use for a solve. Direct maps to a direct
// factorization (SciPP spsolve today), CG to the SPD iterative path (SciPP cg).
// (Same caveats as SciPP#10 apply — see numerics/linear_static.hpp.)
enum class SolverKind {
  Direct,  // scipp::sparse::spsolve (dense LU today — see caveat)
  CG,      // scipp::sparse::cg (unpreconditioned — see caveat)
};

// Which compute backend to run on. CPU is the reference implementation and the
// only one implemented today; the others are reserved for future acceleration and
// never gate a correct build or run.
enum class BackendKind {
  CPU,      // NumPP/SciPP reference path (default, always available)
  CUDA,     // future — not implemented
  OpenCL,   // future — not implemented
  Metal,    // future — not implemented
};

// Abstract compute backend. Backends are stateless value services obtained from
// select_backend(); the interface is intentionally minimal for Phase 1.
class ComputeBackend {
 public:
  virtual ~ComputeBackend() = default;

  // Which backend this instance is.
  virtual BackendKind kind() const = 0;

  // Solve the sparse SPD system A x = rhs, where A is given as COO triplets
  // (rows[k], cols[k], vals[k]); duplicate (row,col) entries are summed. n is the
  // system dimension (rhs.size() == n). Returns x (size n). Uses the requested
  // SolverKind.
  virtual std::vector<Real> solve_sparse(const std::vector<Index>& rows,
                                         const std::vector<Index>& cols,
                                         const std::vector<Real>& vals,
                                         Index n,
                                         const std::vector<Real>& rhs,
                                         SolverKind solver) const = 0;
};

// Stable lowercase name of a backend kind ("cpu", "cuda", "opencl", "metal"),
// suitable for scripting/reporting and round-tripping through backend_kind().
std::string backend_name(BackendKind kind);

// Parse a backend name (case-insensitive; "cpu"/"cuda"/"opencl"/"metal") back to
// a BackendKind. Throws std::invalid_argument on an unknown name.
BackendKind backend_kind(const std::string& name);

// Whether a backend is actually implemented (runnable) on this build. Phase 1: only
// CPU. Unimplemented backends are reported but fall back to CPU when selected.
bool backend_available(BackendKind kind);

// The backends that are actually implemented on this build (CPU only in Phase 1).
std::vector<BackendKind> available_backends();

// The CPU reference backend (always available). Defined in cpu_backend.cpp.
const ComputeBackend& cpu_backend();

// Select a backend. Phase 1 always returns the CPU reference backend regardless
// of the request: unimplemented backends silently fall back to CPU so that the
// absence of a GPU toolkit never breaks the build or run. `requested` documents
// the intended target for callers/tests.
const ComputeBackend& select_backend(BackendKind requested = BackendKind::CPU);

}  // namespace cxpp::compute
