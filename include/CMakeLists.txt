# ==== Define Target func ====
# Assumes these variables are defined:
# TODO

# We shouldn't have written any cpp files by hand so globbing here is fine
file(GLOB func_src ${func_src_dir}/*.cpp)

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

# explicitly declare which headers we're distributing
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
table_types/ArmadilloInterpolationTable.hpp;\
table_types/ConstantTaylorTable.hpp;\
table_types/CubicHermiteTable.hpp;\
table_types/CubicInterpolationTable.hpp;\
table_types/LinearTaylorTable.hpp;\
table_types/QuadraticTaylorTable.hpp;\
table_types/CubicTaylorTable.hpp;\
table_types/LinearInterpolationTable.hpp;\
table_types/LinearRawInterpolationTable.hpp;\
table_types/PadeTable.hpp;\
table_types/QuadraticInterpolationTable.hpp;\
")

target_include_directories(func ${FUNC_PUBLIC}
  $<INSTALL_INTERFACE:.> $<INSTALL_INTERFACE:table_types>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/table_types>
  )

target_link_libraries(func ${FUNC_PUBLIC} ${ARMADILLO_LIBRARIES} ${QUADMATH_LIBRARIES} Boost::boost OpenMP::OpenMP_CXX)

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


# ==== Install Target func ====
include(GNUInstallDirs) # defines CMAKE_INSTALL_XXX
install(TARGETS func
  EXPORT func-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/func
)
set_target_properties(func PROPERTIES EXPORT_NAME FunC)

# make every file in this directory accessable after installing FunC
#install(DIRECTORY ./ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/func)

set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/FunC)
install(EXPORT func-targets
  FILE
    FunCTargets.cmake
  NAMESPACE
    FunC::
  DESTINATION
    ${INSTALL_CONFIGDIR}
)

# Create a ConfigVersion.cmake file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/func-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(${CMAKE_SOURCE_DIR}/cmake/func-config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/func-config.cmake
  INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

install(FILES
  ${CMAKE_SOURCE_DIR}/cmake/FindQuadmath.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/func-config.cmake
  DESTINATION ${INSTALL_CONFIGDIR}
)

# Exporting from the build tree
configure_file(${CMAKE_SOURCE_DIR}/cmake/FindQuadmath.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/FindQuadmath.cmake
  COPYONLY)

export(EXPORT func-targets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/FunCTargets.cmake
  NAMESPACE FunC::)

# Register package in the User Package Registry (TODO is that helpful?)
export(PACKAGE FunC)
