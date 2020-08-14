# FunC library
add_library(func SHARED empty.cpp)

# Implementations are separated out into their own library
list(APPEND func_include_dirs ./ table_types)

if(USE_ARMADILLO AND ARMADILLO_FOUND)
  # link with Armadillo add all the table types that use Armadillo
  list(APPEND func_include_dirs ${ARMADILLO_INCLUDE_DIR})
  list(APPEND func_link_libs ${ARMADILLO_LIBRARIES})
  target_compile_definitions(func PUBLIC -DUSE_ARMADILLO)
endif()

# add the func lib with the lists we built above
target_link_libraries(func PUBLIC ${func_link_libs})
target_include_directories(func PUBLIC ${func_include_dirs})

target_link_libraries(func PUBLIC ${QUADMATH_LIBRARIES} Boost::boost)
if(USE_QUADMATH AND QUADMATH_FOUND)
  target_compile_definitions(func PUBLIC -DUSE_QUADMATH)
endif()

if(FUNC_RECORD)
  target_compile_definitions(func PUBLIC -DFUNC_RECORD)
endif()

set_target_properties(func PROPERTIES PUBLIC_HEADER "\
func.hpp;func_impls.hpp;\
FunctionContainer.hpp;EvaluationImplementation.hpp;\
DirectEvaluation.hpp;ArgumentRecord.hpp;\
ImplementationComparator.hpp;RngInterface.hpp;StdRng.hpp;Timer.hpp;\
json.hpp;\
TransferFunctionInterface.hpp;TransferFunctionSinh.hpp;\
UniformLookupTableGenerator.hpp;\
table_types/TableIncludes.hpp;\
table_types/UniformLookupTable.hpp;\
table_types/UniformArmadilloPrecomputedInterpolationTable.hpp;\
table_types/UniformConstantTaylorTable.hpp;\
table_types/UniformCubicHermiteTable.hpp;\
table_types/UniformCubicPrecomputedInterpolationTable.hpp;\
table_types/UniformCubicTaylorTable.hpp;\
table_types/UniformLinearInterpolationTable.hpp;\
table_types/UniformLinearPrecomputedInterpolationTable.hpp;\
table_types/UniformLinearTaylorTable.hpp;\
table_types/UniformPadeTable.hpp;\
table_types/UniformQuadraticPrecomputedInterpolationTable.hpp;\
table_types/UniformQuadraticTaylorTable.hpp;\
table_types/CompositeLookupTable.hpp;\
table_types/FailureProofTable.hpp;\
table_types/NonUniformLookupTable.hpp;\
table_types/NonUniformLinearInterpolationTable.hpp;\
table_types/NonUniformPseudoLinearInterpolationTable.hpp;\
table_types/NonUniformCubicPrecomputedInterpolationTable.hpp;\
table_types/NonUniformPseudoCubicPrecomputedInterpolationTable.hpp;\
")

install(TARGETS func
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    PUBLIC_HEADER DESTINATION include/func
)
