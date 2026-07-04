# calculixpp_core — dependency-free solver core (domain model + FE element math +
# assembly to COO triplets). No external deps so it builds and tests everywhere,
# including mobile toolchains. The SciPP/NumPP sparse-solve layer links on top.
add_library(calculixpp_core
  core/amplitude.cpp
  core/element.cpp
  core/mesh.cpp
  core/model.cpp
  core/constraints.cpp
  fem/element.cpp
  fem/constraint_transform.cpp
  fem/plasticity.cpp
  fem/hyperelastic.cpp
  fem/user_material.cpp
  fem/assembly.cpp
  fem/stress.cpp
  fem/loads.cpp
  fem/thermal.cpp
  fem/cavity_radiation.cpp
  fem/contact_search.cpp
  fem/contact.cpp
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
    numerics/linear_static.cpp
    numerics/nonlinear_static.cpp
    numerics/heat_transfer.cpp
    numerics/multistep.cpp
    numerics/eigensolution.cpp
    numerics/buckling.cpp
    numerics/modal_dynamics.cpp
    numerics/direct_dynamics.cpp
    numerics/substructure.cpp
    numerics/submodel.cpp
    numerics/high_cycle_fatigue.cpp)
  target_link_libraries(calculixpp_numerics PUBLIC calculixpp_core scipp::scipp)
  target_link_libraries(calculixpp_numerics PRIVATE calculixpp_warnings)
  target_compile_features(calculixpp_numerics PUBLIC cxx_std_20)
  add_library(calculixpp::numerics ALIAS calculixpp_numerics)
endif()
