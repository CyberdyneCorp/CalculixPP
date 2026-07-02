# calculixpp_core — dependency-free solver core (domain model + FE element math +
# assembly to COO triplets). No external deps so it builds and tests everywhere,
# including mobile toolchains. The SciPP/NumPP sparse-solve layer links on top.
add_library(calculixpp_core
  core/element.cpp
  core/mesh.cpp
  core/model.cpp
  fem/element.cpp
  fem/assembly.cpp)

target_include_directories(calculixpp_core
  PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_compile_features(calculixpp_core PUBLIC cxx_std_20)
target_link_libraries(calculixpp_core PRIVATE calculixpp_warnings)

add_library(calculixpp::core ALIAS calculixpp_core)
