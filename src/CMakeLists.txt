# calculixpp_core — dependency-free solver core (domain model + FE element math +
# assembly to COO triplets). No external deps so it builds and tests everywhere,
# including mobile toolchains. The SciPP/NumPP sparse-solve layer links on top.
add_library(calculixpp_core
  core/element.cpp
  core/mesh.cpp
  core/model.cpp
  fem/element.cpp
  fem/assembly.cpp
  fem/stress.cpp
  fem/loads.cpp
  io/inp_parser.cpp
  io/results_writer.cpp
  geometry/cad_mesher.cpp)

target_include_directories(calculixpp_core
  PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_compile_features(calculixpp_core PUBLIC cxx_std_20)
target_link_libraries(calculixpp_core PRIVATE calculixpp_warnings)

add_library(calculixpp::core ALIAS calculixpp_core)

# Numerics layer: sparse solve via SciPP (scipp::sparse) on top of NumPP.
# The CPU reference ComputeBackend lives here since it links SciPP; the backend
# abstraction (compute/backend.hpp) itself is header-only and dependency-free.
if(CALCULIXPP_WITH_SOLVER)
  add_library(calculixpp_numerics
    compute/cpu_backend.cpp
    numerics/linear_static.cpp)
  target_link_libraries(calculixpp_numerics PUBLIC calculixpp_core scipp::scipp)
  target_link_libraries(calculixpp_numerics PRIVATE calculixpp_warnings)
  target_compile_features(calculixpp_numerics PUBLIC cxx_std_20)
  add_library(calculixpp::numerics ALIAS calculixpp_numerics)
endif()
