# build our config file from the given cmake options
file(WRITE config.hpp "#pragma once\n")
if(USE_OPENMP)
  # just a nice option to have for profiling
  # when we start building tables in parallel
  list(APPEND func_link_openmp "OpenMP::OpenMP_CXX")
endif()

if(USE_QUADMATH AND QUADMATH_FOUND)
  file(APPEND config.hpp "#define FUNC_USE_QUADMATH\n")
endif()

if(USE_BOOST_AUTODIFF)  
  file(APPEND config.hpp "#define FUNC_USE_BOOST_AUTODIFF\n")
endif()

# dev option limits our lookup table factory so each
# compilation error only shows up once rather than 4-6 times
if(USE_SMALL_REGISTRY)
  file(APPEND config.hpp "#define FUNC_USE_SMALL_REGISTRY\n")
endif()

# fill our lookup table factory with constructors
list(APPEND func_src
  table_types/StandardRegistrarDefinitions.cpp
)

# Implementations are separated out into their own library
list(APPEND func_include_dirs ./ table_types)

if(USE_ARMADILLO AND ARMADILLO_FOUND)
  # link with Armadillo add all the table types that use Armadillo
  list(APPEND func_include_dirs ${ARMADILLO_INCLUDE_DIR})
  list(APPEND func_link_arma ${ARMADILLO_LIBRARIES})

  # set our config file to use armadillo tables and 
  # fill our lookup table factory with armadillo exclusive table constructors
  file(APPEND config.hpp "#define FUNC_USE_ARMADILLO\n")
  list(APPEND func_src
    table_types/ArmadilloRegistrarDefinitions.cpp
  )
endif()

# add the func lib with the lists we built above
add_library(func SHARED ${func_src})
target_include_directories(func PUBLIC ${func_include_dirs})
target_link_libraries(func PUBLIC ${func_link_arma} ${QUADMATH_LIBRARIES} Boost::boost ${func_link_openmp})

set_target_properties(func PROPERTIES PUBLIC_HEADER "\
config.hpp;func.hpp;\
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
