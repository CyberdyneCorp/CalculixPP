# Phase-1 unit tests. Dependency-free harness (check.hpp) for the first green build;
# GoogleTest integration lands with the wider Phase-1 test task.
foreach(t IN ITEMS test_element test_assembly test_geometry test_step_changes test_contact_search)
  add_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE calculixpp_core)
  target_include_directories(${t} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(NAME ${t} COMMAND ${t})
endforeach()

if(CALCULIXPP_WITH_SOLVER AND TARGET calculixpp_numerics)
  foreach(t IN ITEMS test_solve test_parser test_results test_pressure test_backend test_nonlinear test_loads test_hex test_connector test_plasticity test_material_models test_constraints test_thermal test_thermal_transient test_thermal_stress test_model_change test_cavity_radiation test_multistep test_contact test_thermal_contact test_model_change_contact test_eigen test_buckling test_modal_dynamics test_direct_dynamics test_substructure test_complex_frequency test_submodel test_high_cycle_fatigue)
    add_executable(${t} ${t}.cpp)
    target_link_libraries(${t} PRIVATE calculixpp_numerics)
    target_include_directories(${t} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    add_test(NAME ${t} COMMAND ${t})
  endforeach()
endif()
