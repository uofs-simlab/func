# ==== Define Target func ====
# Assumes these variables are defined:
# TODO

# We shouldn't have written any cpp files by hand so globbing here is fine
file(GLOB func_src *.cpp)

# add the func lib with the lists we built above
if("${func_src}" STREQUAL "")
  message(STATUS "Generating a header-only version of FunC")
  add_library(func INTERFACE)
  set(FUNC_PUBLIC INTERFACE)
  set(FUNC_PRIVATE INTERFACE)
else()
  message(STATUS "Compiling FunC with types ${TYPE_PAIRS}")
  add_library(func SHARED ${func_src})
  set(FUNC_PUBLIC PUBLIC)
  set(FUNC_PRIVATE PRIVATE)
endif()

# TODO I think this needs some fancy generator pathname in case we use make/ make install ???????
target_include_directories(func ${FUNC_PUBLIC} ./ table_types ${ARMADILLO_INCLUDE_DIR}) # TODO should the arma include be there?
target_link_libraries(func ${FUNC_PUBLIC} ${ARMADILLO_LIBRARIES} ${QUADMATH_LIBRARIES} Boost::boost OpenMP::OpenMP_CXX)

set_target_properties(func PROPERTIES PUBLIC_HEADER "\
config.hpp;func.hpp;\
FunctionContainer.hpp;EvaluationImplementation.hpp;\
DirectEvaluation.hpp;ArgumentRecord.hpp;\
ImplementationComparator.hpp;RngInterface.hpp;StdRng.hpp;Timer.hpp;\
json.hpp;\
FailureProofTable.hpp;CompositeLookupTable.hpp;\
TransferFunctionSinh.hpp;\
LookupTableFactory.hpp;\
LookupTableGenerator.hpp;\
table_types/TableIncludes.hpp;\
table_types/LookupTable.hpp;\
table_types/MetaTable.hpp;\
table_types/ArmadilloPrecomputedInterpolationTable.hpp;\
table_types/ConstantTaylorTable.hpp;\
table_types/CubicHermiteTable.hpp;\
table_types/CubicPrecomputedInterpolationTable.hpp;\
table_types/LinearTaylorTable.hpp;\
table_types/QuadraticTaylorTable.hpp;\
table_types/CubicTaylorTable.hpp;\
table_types/LinearInterpolationTable.hpp;\
table_types/LinearPrecomputedInterpolationTable.hpp;\
table_types/PadeTable.hpp;\
table_types/QuadraticPrecomputedInterpolationTable.hpp;\
")

# TODO add more required compile features (it's nice to have cmake determine the required standard on its own)
target_compile_features(func ${FUNC_PRIVATE} cxx_range_for)

# Enable compiler warnings, error on warnings
if (MSVC) # Microsoft
    # warning level 4 and all warnings as errors
    target_compile_options(func ${FUNC_PRIVATE} /W4 /WX)
else() # GCC and Clang
    target_compile_options(func ${FUNC_PRIVATE} -Wall -Wextra -pedantic -Werror)
    # Disable variadic macros warnings (Needed in FUNC_STR_EXPAND)
    target_compile_options(func ${FUNC_PRIVATE} -Wno-variadic-macros)
endif()

# Enable alignment of `new` operator
# TODO does MSVC support this?
target_compile_options(func ${FUNC_PRIVATE} -faligned-new)

install(TARGETS func
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
    PUBLIC_HEADER DESTINATION include/func
)
