### Examples demonstrating the usage of FunC

## Computing maximum errors of all of the implementation types for varying step size
#add_executable(compute_max_err_of_tables compute_max_err_of_tables.cpp)
#target_link_libraries(compute_max_err_of_tables PUBLIC func)

# TODO would it be good to test find_package(FunC) here????
if(FUNC_USE_BOOST)
  add_executable(experiment experiment.cpp)
  target_link_libraries(experiment PUBLIC func)
  target_compile_options(experiment PRIVATE -O3 -Wall -Wextra -pedantic)
endif()

add_executable(plot_impl plot_impl.cpp)
target_link_libraries(plot_impl PUBLIC func)
target_compile_options(plot_impl PRIVATE -O3 -Wall -Wextra -pedantic)

add_executable(compute_max_err_of_tables compute_max_err_of_tables.cpp)
target_link_libraries(compute_max_err_of_tables PUBLIC func)
target_compile_options(compute_max_err_of_tables PRIVATE -O3 -Wall -Wextra -pedantic)

add_executable(test_interfaces test_interfaces.cpp)
target_link_libraries(test_interfaces PUBLIC func)
target_compile_options(test_interfaces PRIVATE -O3 -Wall -Wextra -pedantic)

#add_executable(best-worst best-worst.cpp)
#target_link_libraries(best-worst PUBLIC func)
#target_compile_options(best-worst PRIVATE -O3 -Wall -Wextra -pedantic)

