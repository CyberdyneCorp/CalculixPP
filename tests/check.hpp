#pragma once
// Minimal dependency-free test harness for the first green build.
// (GoogleTest/CTest integration lands with the broader Phase-1 test work.)
#include <cmath>
#include <cstdio>

namespace cxtest {
inline int g_failures = 0;
}

#define CX_CHECK(cond)                                                     \
  do {                                                                     \
    if (!(cond)) {                                                         \
      std::fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
      ++cxtest::g_failures;                                                \
    }                                                                      \
  } while (0)

#define CX_NEAR(a, b, tol)                                                    \
  do {                                                                        \
    const double cx_a = (a), cx_b = (b), cx_t = (tol);                        \
    if (std::fabs(cx_a - cx_b) > cx_t) {                                      \
      std::fprintf(stderr, "FAIL %s:%d: |%.6e - %.6e| = %.3e > %.1e (%s~%s)\n", \
                   __FILE__, __LINE__, cx_a, cx_b, std::fabs(cx_a - cx_b),    \
                   cx_t, #a, #b);                                             \
      ++cxtest::g_failures;                                                   \
    }                                                                         \
  } while (0)

#define CX_MAIN_RETURN() return cxtest::g_failures == 0 ? 0 : 1
