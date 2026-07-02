# Android (NDK) cross-compile toolchain — CONFIGURE-ONLY smoke check.
#
# Purpose: let CMake configure the dependency-free `calculixpp_core` target for
# Android. It is NOT a full build toolchain for the solver: NumPP/SciPP are not
# cross-compiled here, so configure with -DCALCULIXPP_WITH_SOLVER=OFF.
#
# This file defers entirely to the NDK's own `android.toolchain.cmake` for the
# real compiler/sysroot setup, so it tracks whatever NDK the runner provides.
#
# Required: ANDROID_NDK (or ANDROID_NDK_HOME / ANDROID_NDK_ROOT) pointing at an
# installed NDK. Optional cache overrides (standard NDK variables):
#   -DANDROID_ABI=arm64-v8a            (default below)
#   -DANDROID_PLATFORM=android-24      (default below)
#
# Usage:
#   cmake -S . -B build-android \
#     -DCMAKE_TOOLCHAIN_FILE=cmake/toolchains/android.cmake \
#     -DCALCULIXPP_WITH_SOLVER=OFF -DCALCULIXPP_BUILD_TESTS=OFF

# Locate the NDK from the common environment variables if not passed explicitly.
if(NOT DEFINED ANDROID_NDK)
  if(DEFINED ENV{ANDROID_NDK})
    set(ANDROID_NDK "$ENV{ANDROID_NDK}")
  elseif(DEFINED ENV{ANDROID_NDK_HOME})
    set(ANDROID_NDK "$ENV{ANDROID_NDK_HOME}")
  elseif(DEFINED ENV{ANDROID_NDK_ROOT})
    set(ANDROID_NDK "$ENV{ANDROID_NDK_ROOT}")
  endif()
endif()

if(NOT ANDROID_NDK OR NOT EXISTS "${ANDROID_NDK}")
  message(FATAL_ERROR
    "Android NDK not found. Set ANDROID_NDK (or the ANDROID_NDK_HOME / "
    "ANDROID_NDK_ROOT environment variable) to an installed NDK.")
endif()

# Standard NDK knobs — sensible defaults, overridable from the command line.
if(NOT DEFINED ANDROID_ABI)
  set(ANDROID_ABI "arm64-v8a")
endif()
if(NOT DEFINED ANDROID_PLATFORM)
  set(ANDROID_PLATFORM "android-24")
endif()

# Hand off to the NDK's own toolchain, which sets CMAKE_SYSTEM_NAME=Android,
# the clang cross-compiler, sysroot, and CMAKE_ANDROID_* variables.
include("${ANDROID_NDK}/build/cmake/android.toolchain.cmake")
