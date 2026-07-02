# Python bindings (pybind11). Locate pybind11's CMake package via the Python
# interpreter so it works from a pip/conda install without extra flags.
find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

if(NOT pybind11_DIR)
  execute_process(
    COMMAND "${Python3_EXECUTABLE}" -m pybind11 --cmakedir
    OUTPUT_VARIABLE _pybind11_dir OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE _pybind11_rc)
  if(_pybind11_rc EQUAL 0)
    set(pybind11_DIR "${_pybind11_dir}")
  endif()
endif()

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(calculixpp_py bindings.cpp)
set_target_properties(calculixpp_py PROPERTIES OUTPUT_NAME calculixpp)
target_link_libraries(calculixpp_py PRIVATE calculixpp_numerics)

# Regression harness: solve reference decks and diff against CalculiX output.
if(CALCULIXPP_BUILD_TESTS)
  add_test(NAME python_regression
    COMMAND "${Python3_EXECUTABLE}" -m pytest -q
            "${CMAKE_CURRENT_SOURCE_DIR}/tests")
  set_tests_properties(python_regression PROPERTIES
    ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:calculixpp_py>")
endif()
