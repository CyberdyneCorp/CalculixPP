# Phase-1 unit tests. Dependency-free harness (check.hpp) for the first green build;
# GoogleTest integration lands with the wider Phase-1 test task.
foreach(t IN ITEMS test_element test_assembly)
  add_executable(${t} ${t}.cpp)
  target_link_libraries(${t} PRIVATE calculixpp_core)
  target_include_directories(${t} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(NAME ${t} COMMAND ${t})
endforeach()

if(CALCULIXPP_WITH_SOLVER AND TARGET calculixpp_numerics)
  add_executable(test_solve test_solve.cpp)
  target_link_libraries(test_solve PRIVATE calculixpp_numerics)
  target_include_directories(test_solve PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  add_test(NAME test_solve COMMAND test_solve)
endif()
