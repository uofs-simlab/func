
# Implementations are separated out into their own library
# Create lists for func_impls source files, libs to link with, and include
# dir in case the user doesn't have Armadillo installed
list(APPEND func_impls_src
    EvaluationImplementation.cpp
    DirectEvaluation.cpp
    ArgumentRecord.cpp
    table_types/UniformConstantTaylorTable.cpp
    table_types/UniformCubicHermiteTable.cpp
    table_types/UniformCubicPrecomputedInterpolationTable.cpp
    table_types/UniformCubicTaylorTable.cpp
    table_types/UniformFailureProofCubicPITable.cpp
    table_types/UniformFailureProofLinearPITable.cpp
    table_types/UniformLinearInterpolationTable.cpp
    table_types/UniformLinearPrecomputedInterpolationTable.cpp
    table_types/UniformLinearTaylorTable.cpp
    table_types/UniformLookupTable.cpp
    table_types/UniformQuadraticPrecomputedInterpolationTable.cpp
    table_types/UniformQuadraticTaylorTable.cpp
)
list(APPEND func_impls_include_dirs ./ table_types)

if(USE_ARMADILLO AND ARMADILLO_FOUND)
  # link with Armadillo add all the table types that use Armadillo
  list(APPEND func_impls_include_dirs ${ARMADILLO_INCLUDE_DIR})
  list(APPEND func_impls_link_libs ${ARMADILLO_LIBRARIES})

  list(APPEND func_impls_src
      table_types/UniformArmadilloPrecomputedInterpolationTable.cpp
      table_types/UniformPadeTable.cpp
  )
endif()

# add the func_impls lib with the lists we built above
add_library(func_impls SHARED ${func_impls_src})
target_link_libraries(func_impls PUBLIC ${func_impls_link_libs})
target_include_directories(func_impls PUBLIC ${func_impls_include_dirs})

# FunC library
add_library(func SHARED
    ImplementationComparator.cpp
    Timer.cpp
    UniformLookupTableGenerator.cpp
)

target_link_libraries(func PUBLIC ${QUADMATH_LIBRARIES} Boost::boost func_impls)
if(USE_QUADMATH AND QUADMATH_FOUND)
  target_compile_definitions(func PUBLIC -DUSE_QUADMATH)
endif()

if(FUNC_RECORD)
  target_compile_definitions(func PUBLIC -DFUNC_RECORD)
  target_compile_definitions(func_impls PUBLIC -DFUNC_RECORD)
endif()

if(NDEBUG)
  target_compile_definitions(func_impls PUBLIC -NDEBUG)
endif()

# TODO: Create a simplified public header? BTW every header is always being made
set_target_properties(func PROPERTIES PUBLIC_HEADER "func.hpp;func_impls.hpp;EvaluationFunctor.hpp;EvaluationImplementation.hpp;DirectEvaluation.hpp;ImplementationComparator.hpp;RngInterface.hpp;StdRng.hpp;Timer.hpp;UniformLookupTableGenerator.hpp;table_types/UniformArmadilloPrecomputedInterpolationTable.hpp;table_types/UniformConstantTaylorTable.hpp;table_types/UniformCubicHermiteTable.hpp;table_types/UniformCubicPrecomputedInterpolationTable.hpp;table_types/UniformCubicTaylorTable.hpp;table_types/UniformFailureProofCubicPITable.hpp;table_types/UniformFailureProofLinearPITable.hpp;table_types/UniformLinearInterpolationTable.hpp;table_types/UniformLinearPrecomputedInterpolationTable.hpp;table_types/UniformLinearTaylorTable.hpp;table_types/UniformLookupTable.hpp;table_types/UniformPadeTable.hpp;table_types/UniformQuadraticPrecomputedInterpolationTable.hpp;table_types/UniformQuadraticTaylorTable.hpp;table_types/UniformTables.hpp")

install(TARGETS func func_impls
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    PUBLIC_HEADER DESTINATION include/func
)
