# Numerics dependencies for the solver layer.
#
# NumPP ships a CMake config package -> find_package(NumPP).
# SciPP has no install/export -> consumed via add_subdirectory (it internally does
# find_package(NumPP CONFIG REQUIRED), resolved from the same prefix). SciPP provides
# the sparse module; NumPP provides ndarray + dense linear algebra.
#
# Run scripts/bootstrap_deps.sh first to build+install NumPP into .deps/install.

set(CALCULIXPP_SCIPP_DIR "${CMAKE_SOURCE_DIR}/../SciPP"
    CACHE PATH "SciPP source checkout")

find_package(NumPP CONFIG REQUIRED)

if(NOT TARGET scipp::scipp)
  if(NOT EXISTS "${CALCULIXPP_SCIPP_DIR}/CMakeLists.txt")
    message(FATAL_ERROR
      "SciPP not found. Set -DCALCULIXPP_SCIPP_DIR=<path-to-SciPP> "
      "(looked in '${CALCULIXPP_SCIPP_DIR}').")
  endif()
  # We only need the SciPP library, not its own test suite.
  set(SCIPP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  add_subdirectory("${CALCULIXPP_SCIPP_DIR}" "${CMAKE_BINARY_DIR}/_scipp")
endif()

message(STATUS "CalculiX++ solver deps: NumPP ${NumPP_VERSION}; SciPP from ${CALCULIXPP_SCIPP_DIR}")
