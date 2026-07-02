#include "calculixpp/core/amplitude.hpp"

#include <algorithm>
#include <cmath>

namespace cxpp {

Real Amplitude::value_at(Real t) const {
  if (points.empty()) return 1.0;
  const Real t0 = points.front().first;

  // Wrap a periodic amplitude into [t0, t0 + period).
  if (period > 0.0) {
    const Real rel = t - t0;
    const Real wrapped = rel - std::floor(rel / period) * period;
    t = t0 + wrapped;
  }

  if (t <= t0) return points.front().second;
  const Real tN = points.back().first;
  if (t >= tN) return points.back().second;

  // Find the bracketing interval [i, i+1] with time_i <= t < time_{i+1}.
  std::size_t i = 0;
  while (i + 1 < points.size() && points[i + 1].first <= t) ++i;

  if (definition == Definition::Step) return points[i].second;

  const auto& [ta, va] = points[i];
  const auto& [tb, vb] = points[i + 1];
  const Real span = tb - ta;
  if (span <= 0.0) return vb;  // degenerate interval: take the later value
  return va + (vb - va) * (t - ta) / span;
}

}  // namespace cxpp
