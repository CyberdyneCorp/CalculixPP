#include "calculixpp/compute/backend.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <stdexcept>

#include "numpp/core/creation.hpp"
#include "numpp/core/dtype.hpp"
#include "scipp/sparse/sparse.hpp"

// CPU reference backend (spec: compute-backend, linear-algebra-and-solvers).
// Wraps the SciPP sparse path: COO triplets -> CsrMatrix::from_coo -> spsolve
// (sparse Cholesky for SPD / LU otherwise, SciPP v1.2.0) or IC0-preconditioned
// CG. Direct is the default; CG uses cg_report + IC0 and checks convergence so a
// non-converged result raises rather than returning a silently-wrong vector.
namespace cxpp::compute {
namespace {

class CpuBackend final : public ComputeBackend {
 public:
  BackendKind kind() const override { return BackendKind::CPU; }

  std::vector<Real> solve_sparse(const std::vector<Index>& rows,
                                 const std::vector<Index>& cols,
                                 const std::vector<Real>& vals,
                                 Index n,
                                 const std::vector<Real>& rhs,
                                 SolverKind solver) const override {
    const std::int64_t dim = static_cast<std::int64_t>(n);
    const std::int64_t nnz = static_cast<std::int64_t>(vals.size());

    scipp::sparse::CooMatrix coo;
    coo.rows = dim;
    coo.cols = dim;
    coo.data = numpp::zeros({nnz}, numpp::kFloat64);
    coo.row = numpp::zeros({nnz}, numpp::kInt64);
    coo.col = numpp::zeros({nnz}, numpp::kInt64);
    for (std::int64_t i = 0; i < nnz; ++i) {
      const auto k = static_cast<std::size_t>(i);
      coo.data.set_item<double>({i}, vals[k]);
      coo.row.set_item<std::int64_t>({i}, rows[k]);
      coo.col.set_item<std::int64_t>({i}, cols[k]);
    }

    const scipp::sparse::CsrMatrix A = scipp::sparse::CsrMatrix::from_coo(coo);

    numpp::ndarray b = numpp::zeros({dim}, numpp::kFloat64);
    for (std::int64_t i = 0; i < dim; ++i)
      b.set_item<double>({i}, rhs[static_cast<std::size_t>(i)]);

    numpp::ndarray x;
    if (solver == SolverKind::CG) {
      // IC0-preconditioned CG for the SPD FE system, with a convergence check.
      const scipp::sparse::IterationResult res =
          scipp::sparse::cg_report(A, b, scipp::sparse::Preconditioner::IC0);
      if (!res.converged)
        throw std::runtime_error(
            "CG did not converge (relative residual " +
            std::to_string(res.final_residual) + " after " +
            std::to_string(res.iterations) +
            " iterations); use SOLVER=default (direct)");
      x = res.x;
    } else {
      x = scipp::sparse::spsolve(A, b);  // sparse Cholesky (SPD) / LU
    }

    std::vector<Real> u(static_cast<std::size_t>(dim));
    for (std::int64_t i = 0; i < dim; ++i)
      u[static_cast<std::size_t>(i)] = x.item<double>({i});
    return u;
  }
};

}  // namespace

std::string backend_name(BackendKind kind) {
  switch (kind) {
    case BackendKind::CPU: return "cpu";
    case BackendKind::CUDA: return "cuda";
    case BackendKind::OpenCL: return "opencl";
    case BackendKind::Metal: return "metal";
  }
  return "cpu";
}

BackendKind backend_kind(const std::string& name) {
  std::string lower = name;
  std::transform(lower.begin(), lower.end(), lower.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  if (lower == "cpu") return BackendKind::CPU;
  if (lower == "cuda") return BackendKind::CUDA;
  if (lower == "opencl") return BackendKind::OpenCL;
  if (lower == "metal") return BackendKind::Metal;
  throw std::invalid_argument("unknown compute backend: '" + name +
                              "' (expected cpu/cuda/opencl/metal)");
}

bool backend_available(BackendKind kind) {
  // Phase 1: only the CPU reference backend is implemented.
  return kind == BackendKind::CPU;
}

std::vector<BackendKind> available_backends() {
  return {BackendKind::CPU};
}

const ComputeBackend& cpu_backend() {
  static const CpuBackend backend;
  return backend;
}

const ComputeBackend& select_backend(BackendKind /*requested*/) {
  // Only the CPU backend is implemented in Phase 1. Any request falls back to CPU
  // so that a missing GPU toolkit never breaks the build or run.
  return cpu_backend();
}

}  // namespace cxpp::compute
