# FunC libs
add_library(func INTERFACE)

# Implementations are separated out into their own library
add_library(func_impls INTERFACE)

list(APPEND func_impls_headers
    ${CMAKE_CURRENT_SOURCE_DIR}/func_impls.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FunctionContainer.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/EvaluationImplementation.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DirectEvaluation.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/table_types/TableIncludes.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/table_types/UniformLookupTable.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/table_types/UniformConstantTaylorTable.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/table_types/UniformLinearTaylorTable.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/table_types/UniformQuadraticTaylorTable.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/table_types/UniformCubicTaylorTable.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/table_types/UniformCubicHermiteTable.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/table_types/UniformLinearInterpolationTable.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/table_types/UniformLinearPrecomputedInterpolationTable.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/table_types/UniformCubicPrecomputedInterpolationTable.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/table_types/UniformQuadraticPrecomputedInterpolationTable.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/table_types/FailureProofTable.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/table_types/CompositeLookupTable.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TransferFunction.hpp
)

list(APPEND func_impls_include_dirs ./ table_types)

if(USE_ARMADILLO AND ARMADILLO_FOUND)
  # link with Armadillo add all the table types that use Armadillo
  list(APPEND func_impls_include_dirs ${ARMADILLO_INCLUDE_DIR})
  list(APPEND func_impls_link_libs ${ARMADILLO_LIBRARIES})
  target_compile_definitions(func_impls INTERFACE -DUSE_ARMADILLO)
  list(APPEND func_impls_headers 
    ${CMAKE_CURRENT_SOURCE_DIR}/table_types/UniformArmadilloPrecomputedInterpolationTable.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/table_types/UniformPadeTable.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TransferFunctionSinh.hpp
  )
endif()

target_link_libraries(func_impls INTERFACE ${func_impls_link_libs})
target_include_directories(func_impls INTERFACE ${func_impls_include_dirs})
target_sources(func_impls INTERFACE ${func_impls_headers})

#NonUniformLookupTable.hpp;\
#NonUniformLinearInterpolationTable.hpp;\

target_link_libraries(func INTERFACE ${QUADMATH_LIBRARIES} Boost::boost func_impls)
if(USE_QUADMATH AND QUADMATH_FOUND)
  target_compile_definitions(func INTERFACE -DUSE_QUADMATH)
endif()

if(FUNC_RECORD)
  target_compile_definitions(func INTERFACE -DFUNC_RECORD)
  target_compile_definitions(func_impls INTERFACE -DFUNC_RECORD)
endif()

target_sources(func INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/func.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/ImplementationComparator.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/RngInterface.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/StdRng.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Timer.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/UniformLookupTableGenerator.hpp
)

install(TARGETS func func_impls
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    PUBLIC_HEADER DESTINATION include/func
)
