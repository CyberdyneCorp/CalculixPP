cmake_minimum_required(VERSION 3.24)
project(CalculixPP VERSION 0.1.0 LANGUAGES CXX)

# ---------------------------------------------------------------------------
# CalculiX++ — modern C++20 port of the CalculiX finite element solver.
# See openspec/ for the specification and phased roadmap.
# ---------------------------------------------------------------------------

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "" FORCE)
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(CALCULIXPP_BUILD_TESTS  "Build the test suite"                 ON)
option(CALCULIXPP_WITH_SOLVER  "Build the SciPP/NumPP sparse solver layer" ON)
option(CALCULIXPP_WARNINGS_AS_ERRORS "Treat warnings as errors on the core" OFF)

# Portability guardrail (spec: build-and-tooling / compute-backend):
# a correct build requires NO GPU toolkit. GPU backends are additive, later.

add_library(calculixpp_warnings INTERFACE)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(calculixpp_warnings INTERFACE
    -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion)
  if(CALCULIXPP_WARNINGS_AS_ERRORS)
    target_compile_options(calculixpp_warnings INTERFACE -Werror)
  endif()
endif()

add_subdirectory(src)

if(CALCULIXPP_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

message(STATUS "CalculiX++ ${PROJECT_VERSION}  tests=${CALCULIXPP_BUILD_TESTS}  solver=${CALCULIXPP_WITH_SOLVER}")
