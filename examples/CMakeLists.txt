### Examples demonstrating the usage of FunC

## Computing maximum errors of all of the implementation types for varying step size
#add_executable(compute_max_err_of_tables compute_max_err_of_tables.cpp)
#target_link_libraries(compute_max_err_of_tables PUBLIC func)

# TODO is it possible to test find_package(FunC) here????
if(FUNC_USE_BOOST)
  add_executable(experiment experiment.cpp)
  target_link_libraries(experiment PUBLIC func)
  target_compile_options(experiment PRIVATE -O3 -Wall -Wextra -pedantic)
endif()

add_executable(plot_impl plot_impl.cpp)
target_link_libraries(plot_impl PUBLIC func)
target_compile_options(plot_impl PRIVATE -O3 -Wall -Wextra -pedantic)

#add_executable(compute_max_err_of_tables compute_max_err_of_tables.cpp)
#target_link_libraries(compute_max_err_of_tables PUBLIC func)
#target_compile_options(compute_max_err_of_tables PRIVATE -O3 -Wall -Wextra -pedantic)

add_executable(test_interfaces test_interfaces.cpp)
target_link_libraries(test_interfaces PUBLIC func)
target_compile_options(test_interfaces PRIVATE -O3 -Wall -Wextra -pedantic)

add_executable(difftables difftables.cpp)
target_link_libraries(difftables PUBLIC func)
target_compile_options(difftables PRIVATE -O3 -Wall -Wextra -pedantic)

if(FUNC_NO_CXX17)
else()
  add_executable(2D_lut 2D_lut.cpp)
  target_link_libraries(2D_lut PUBLIC func)
  target_compile_options(2D_lut PRIVATE -O3 -Wall -Wextra -pedantic)
endif()


add_test(NAME BuildTables COMMAND experiment 0.001 30.0 1e-2 10 1000000 19)
add_test(NAME TestUtils COMMAND test_interfaces)
# TODO Use CMake to check for the user's installation of python
#add_test(NAME PlotIt COMMAND ./examples/plot_impl "UniformChebyInterpTable<1>" 0.1 30 5 > foo.txt && python3.11 plot.py foo.txt)

#add_executable(best-worst best-worst.cpp)
#target_link_libraries(best-worst PUBLIC func)
#target_compile_options(best-worst PRIVATE -O3 -Wall -Wextra -pedantic)

