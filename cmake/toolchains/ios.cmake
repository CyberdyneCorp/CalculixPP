# iOS cross-compile toolchain — CONFIGURE-ONLY smoke check.
#
# Purpose: let CMake configure the dependency-free `calculixpp_core` target for
# iOS. It is NOT a full build toolchain for the solver: NumPP/SciPP are not
# cross-compiled here, so configure with -DCALCULIXPP_WITH_SOLVER=OFF.
#
# Relies on CMake's built-in iOS support (CMAKE_SYSTEM_NAME=iOS), which drives
# the Xcode toolchain / clang with the iphoneos SDK. Requires a macOS host with
# Xcode + command-line tools installed.
#
# Optional cache overrides:
#   -DCMAKE_OSX_SYSROOT=iphoneos              (default; use iphonesimulator for sim)
#   -DCMAKE_OSX_DEPLOYMENT_TARGET=13.0        (default below)
#   -DCMAKE_OSX_ARCHITECTURES=arm64           (default below)
#
# Usage (on macOS):
#   cmake -S . -B build-ios \
#     -DCMAKE_TOOLCHAIN_FILE=cmake/toolchains/ios.cmake \
#     -DCALCULIXPP_WITH_SOLVER=OFF -DCALCULIXPP_BUILD_TESTS=OFF

# CMAKE_SYSTEM_NAME being set here marks this as a cross-compile and enables
# CMake's Platform/iOS support module.
set(CMAKE_SYSTEM_NAME iOS)

if(NOT DEFINED CMAKE_OSX_SYSROOT)
  set(CMAKE_OSX_SYSROOT iphoneos CACHE STRING "iOS SDK sysroot")
endif()
if(NOT DEFINED CMAKE_OSX_DEPLOYMENT_TARGET)
  set(CMAKE_OSX_DEPLOYMENT_TARGET "13.0" CACHE STRING "Minimum iOS version")
endif()
if(NOT DEFINED CMAKE_OSX_ARCHITECTURES)
  set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "iOS target architectures")
endif()

# Build static libraries only for the configure smoke check; iOS forbids the
# unversioned shared libraries CMake would otherwise try to produce.
set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)
