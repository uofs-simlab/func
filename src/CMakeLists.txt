
# Implementations are separated out into their own library
add_library(func_impls SHARED
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
target_include_directories(func_impls PUBLIC ./ table_types)

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

if(USE_ARMADILLO AND ARMADILLO_FOUND)
  # Compile all the table types that use Armadillo
  target_include_directories(func_impls PUBLIC ${ARMADILLO_INCLUDE_DIR})
  target_link_libraries(func_impls PUBLIC ${ARMADILLO_LIBRARIES})
  add_library(APPEND func_impls
    table_types/UniformArmadilloPrecomputedInterpolationTable.cpp
    table_types/UniformPadeTable.cpp
  )
endif()

if(FUNC_RECORD)
  target_compile_definitions(func PUBLIC -DFUNC_RECORD)
  target_compile_definitions(func_impls PUBLIC -DFUNC_RECORD)
endif()

if(NDEBUG)
  target_compile_definitions(func_impls PUBLIC -NDEBUG)
endif()

# TODO: Create a simplified public header?
set_target_properties(func PROPERTIES PUBLIC_HEADER "func.hpp;func_impls.hpp;EvaluationFunctor.hpp;EvaluationImplementation.hpp;DirectEvaluation.hpp;ImplementationComparator.hpp;RngInterface.hpp;StdRng.hpp;Timer.hpp;UniformLookupTableGenerator.hpp;table_types/UniformArmadilloPrecomputedInterpolationTable.hpp;table_types/UniformConstantTaylorTable.hpp;table_types/UniformCubicHermiteTable.hpp;table_types/UniformCubicPrecomputedInterpolationTable.hpp;table_types/UniformCubicTaylorTable.hpp;table_types/UniformFailureProofCubicPITable.hpp;table_types/UniformFailureProofCubicPITable.hpp;table_types/UniformLinearInterpolationTable.hpp;table_types/UniformLinearPrecomputedInterpolationTable.hpp;table_types/UniformLinearTaylorTable.hpp;table_types/UniformLookupTable.hpp;table_types/UniformPadeTable.hpp;table_types/UniformQuadraticPrecomputedInterpolationTable.hpp;table_types/UniformQuadraticTaylorTable.hpp;table_types/UniformTables.hpp")

install(TARGETS func func_impls
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    PUBLIC_HEADER DESTINATION include/func
)
