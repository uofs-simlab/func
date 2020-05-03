### Examples demonstrating the usage of FunC

## Computing maximum errors of all of the implementation types for varying step size
add_executable(compute_max_err_of_tables compute_max_err_of_tables.cpp)
target_link_libraries(compute_max_err_of_tables PUBLIC func)

add_executable(experiment experiment.cpp)
target_link_libraries(experiment PUBLIC func)

add_executable(experiment_best_worst experiment_best_worst.cpp)
target_link_libraries(experiment_best_worst PUBLIC func)

add_executable(generate_tables_at_tol generate_tables_at_tol.cpp)
target_link_libraries(generate_tables_at_tol PUBLIC func)
