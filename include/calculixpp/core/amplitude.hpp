#pragma once
#include <string>
#include <utility>
#include <vector>

#include "calculixpp/core/types.hpp"

namespace cxpp {

// A time-value amplitude curve (*AMPLITUDE). The step time is normalized to the
// incrementation total; loads/BCs that reference an amplitude are scaled by
// amplitude(name, t) instead of the default linear 0->1 ramp (spec:
// nonlinear-solution-control — Amplitude-driven time stepping; ref: amplitudes.f).
struct Amplitude {
  // How the curve behaves between/around the tabulated points.
  enum class Definition {
    Tabular,  // piecewise-linear interpolation between points
    Step,     // piecewise-constant (holds the last point's value)
  };

  std::string name;
  Definition definition{Definition::Tabular};
  // (time, value) pairs, kept in the order given (assumed non-decreasing time).
  std::vector<std::pair<Real, Real>> points;
  // Period for a periodic amplitude (TIME= / period, > 0); 0 means non-periodic.
  Real period{0.0};

  // Value of the curve at time t. Empty curve -> 1.0 (no scaling). Times outside
  // the tabulated range clamp to the first/last value; when `period` > 0 the time
  // is wrapped into [t0, t0+period) first.
  Real value_at(Real t) const;
};

}  // namespace cxpp
