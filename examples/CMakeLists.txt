### Examples demonstrating the usage of FunC

## Computing maximum errors of all of the implementation types for varying step size
#add_executable(compute_max_err_of_tables compute_max_err_of_tables.cpp)
#target_link_libraries(compute_max_err_of_tables PUBLIC func)

# TODO would it be good to test find_package(FunC) here????
if(FUNC_USE_BOOST)
add_executable(experiment experiment.cpp)
target_link_libraries(experiment PUBLIC func)
target_compile_options(experiment PRIVATE -Wall -Wextra -pedantic)
endif()

add_executable(plot_impl plot_impl.cpp)
target_link_libraries(plot_impl PUBLIC func)
target_compile_options(plot_impl PRIVATE -Wall -Wextra -pedantic)

add_executable(experiment_best_worst experiment_best_worst.cpp)
target_link_libraries(experiment_best_worst PUBLIC func)
target_compile_options(experiment_best_worst PRIVATE -Wall -Wextra -pedantic)

add_executable(test_interfaces test_interfaces.cpp)
target_link_libraries(test_interfaces PUBLIC func)
target_compile_options(test_interfaces PRIVATE -Wall -Wextra -pedantic)

#add_executable(generate_tables_at_tol generate_tables_at_tol.cpp)
#target_link_libraries(generate_tables_at_tol PUBLIC func)
