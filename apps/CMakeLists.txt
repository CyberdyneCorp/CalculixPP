# ccxpp CLI runner (needs the solver layer).
add_executable(ccxpp cli/main.cpp)
target_link_libraries(ccxpp PRIVATE calculixpp_numerics)
