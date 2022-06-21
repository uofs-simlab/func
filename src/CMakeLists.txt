# build our config file from the given cmake options
file(WRITE config.hpp "#pragma once\n")
if(USE_BOOST)
  file(APPEND config.hpp "#define FUNC_USE_BOOST\n")
  list(APPEND func_link_boost Boost::boost)
endif()

if(USE_OPENMP)
  # just a nice option to have for profiling
  # when we start building tables in parallel
  list(APPEND func_link_openmp "OpenMP::OpenMP_CXX")
endif()

if(USE_QUADMATH)
  file(APPEND config.hpp "#define FUNC_USE_QUADMATH\n")
  list(APPEND func_link_quadmath ${QUADMATH_LIBRARIES} )
endif()

# dev option limits our lookup table factory so each
# compilation error only shows up once rather than 4-6 times
if(USE_SMALL_REGISTRY)
  file(APPEND config.hpp "#define FUNC_USE_SMALL_REGISTRY\n")
endif()

if(USE_ARMADILLO)
  # link with Armadillo add all the table types that use Armadillo
  list(APPEND func_include_dirs ${ARMADILLO_INCLUDE_DIR})
  list(APPEND func_link_arma ${ARMADILLO_LIBRARIES})

  # set our config file to use armadillo tables and
  # fill our lookup table factory with armadillo exclusive table constructors
  file(APPEND config.hpp "#define FUNC_USE_ARMADILLO\n")
endif()

# Get all of the cpp files in table_types
# TODO: Separate out table types again (parallel builds)
# TODO: Re-add default lib builds for: double,double  ... float,float ... etc
#       - maybe by forcing UniformLookupTableGenerator<TIN,TOUT> builds?
file(GLOB func_src table_types/Registrar*.cpp)
# Add to include directory
# TODO: Re-add func_impls as its own separate library (again, parallel builds)
list(APPEND func_include_dirs ./ table_types)

# add the func lib with the lists we built above
add_library(func SHARED ${func_src})
target_include_directories(func PUBLIC ${func_include_dirs})
target_link_libraries(func PUBLIC ${func_link_arma} ${func_link_quadmath} Boost::boost ${func_link_openmp})

set_target_properties(func PROPERTIES PUBLIC_HEADER "\
config.hpp;func.hpp;\
FunctionContainer.hpp;EvaluationImplementation.hpp;\
DirectEvaluation.hpp;ArgumentRecord.hpp;\
ImplementationComparator.hpp;RngInterface.hpp;StdRng.hpp;Timer.hpp;\
json.hpp;\
FailureProofTable.hpp;CompositeLookupTable.hpp;\
TransferFunctionInterface.hpp;TransferFunctionSinh.hpp;\
LookupTableFactory.hpp;\
LookupTableGenerator.hpp;\
table_types/TableIncludes.hpp;\
table_types/LookupTable.hpp;\
table_types/MetaTable.hpp;\
table_types/ArmadilloPrecomputedInterpolationTable.hpp;\
table_types/ConstantTaylorTable.hpp;\
table_types/CubicHermiteTable.hpp;\
table_types/CubicPrecomputedInterpolationTable.hpp;\
table_types/CubicTaylorTable.hpp;\
table_types/LinearInterpolationTable.hpp;\
table_types/LinearPrecomputedInterpolationTable.hpp;\
table_types/LinearTaylorTable.hpp;\
table_types/PadeTable.hpp;\
table_types/QuadraticPrecomputedInterpolationTable.hpp;\
table_types/QuadraticTaylorTable.hpp;\
")

install(TARGETS func
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    PUBLIC_HEADER DESTINATION include/func
)
