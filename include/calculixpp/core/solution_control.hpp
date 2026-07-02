#pragma once
#include <vector>

#include "calculixpp/core/types.hpp"

namespace cxpp {

// Newton-Raphson convergence controls, parsed from *CONTROLS (spec:
// nonlinear-solution-control — Convergence controls). Defaults follow CalculiX's
// documented values (controlss.f) closely enough for the linear reproduction gate.
struct NonlinearControls {
  // Force-residual convergence tolerance: ||r|| / max(||f_ext||, eps) < force_tol.
  Real force_tol{5e-3};
  // Displacement-correction tolerance: ||du|| / max(||u||, eps) < disp_tol.
  Real disp_tol{1e-3};
  // Maximum Newton iterations per increment before a cutback is requested.
  int max_iterations{16};
  // Floor for the residual denominator to avoid division by zero on a zero load.
  Real eps{1e-30};
};

// Automatic-incrementation parameters, parsed from the *STATIC data line and its
// DIRECT parameter (spec: nonlinear-solution-control — Automatic incrementation and
// cutback). The step spans load factor 0 -> 1 (total_time normalized internally).
struct Incrementation {
  Real initial{1.0};   // initial increment size (fraction of the step)
  Real total{1.0};     // total step time; increments are normalized against this
  Real min{1e-5};      // minimum increment; the step aborts below this after cutback
  Real max{1.0};       // maximum increment size
  bool direct{false};  // DIRECT: fixed user increment, no automatic resizing
  Real cutback{0.25};  // factor applied to the increment on non-convergence
  Real grow{1.5};      // factor applied after an easily-converged increment
  int grow_below{5};   // grow only if the increment converged in <= this many iters
  int max_increments{1000};  // safety cap on the number of increments in a step
};

// Ordered list of *TIME POINTS (as fractions of the step), used to force increments
// to land exactly on the given times (spec: nonlinear-solution-control).
struct TimePoints {
  std::vector<Real> times;  // strictly increasing, in (0, total]
};

}  // namespace cxpp
