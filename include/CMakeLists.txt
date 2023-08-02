# ==== Define Target func ====
# Assumes these variables are defined:
# ARMADILLO_LIBRARIES QUADMATH_LIBRARIES Boost_LIBRARIES OpenMP_CXX_LIBRARIES

# We didn't write any cpp files by hand so it's good to glob here
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
FunctionContainer.hpp;LookupTable.hpp;config.hpp;func.hpp;json.hpp;\
tables/tables.hpp;utils/utils.hpp;\
tables/TransferFunction.hpp;tables/Polynomial.hpp;tables/MetaTable.hpp;\
tables/ChebyInterpTable.hpp;\
tables/CubicHermiteTable.hpp;\
tables/EqSpaceInterpTable.hpp;\
tables/LinearRawInterpTable.hpp;\
tables/PadeTable.hpp;\
tables/TaylorTable.hpp;\
utils/DirectEvaluation.hpp;utils/ArgumentRecord.hpp;\
utils/LookupTableFactory.hpp;utils/LookupTableGenerator.hpp;\
utils/FailureProofTable.hpp;utils/CompositeLookupTable.hpp;\
utils/RngInterface.hpp;utils/StdRng.hpp;\
utils/LookupTableComparator.hpp;\
utils/Timer.hpp;\
")

target_include_directories(func ${FUNC_PUBLIC}
  $<INSTALL_INTERFACE:include>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/tables> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/utils>
)

target_link_libraries(func ${FUNC_PUBLIC} ${ARMADILLO_LIBRARIES} ${Boost_HEADERS_LIBRARY} ${OPENMP_CXX_LIBRARIES})

# Enable compiler warnings, error on warnings
if (MSVC) # Microsoft
    # warning level 4 and all warnings as errors
    target_compile_options(func ${FUNC_PRIVATE} /W4 /WX)
else() # GCC and Clang
    target_compile_options(func ${FUNC_PRIVATE} -Wall -Wextra -pedantic) #-Werror)
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
  #${CMAKE_SOURCE_DIR}/cmake/FindQuadmath.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/func-config.cmake
  DESTINATION ${INSTALL_CONFIGDIR}
)

# Exporting from the build tree
#configure_file(${CMAKE_SOURCE_DIR}/cmake/FindQuadmath.cmake
#  ${CMAKE_CURRENT_BINARY_DIR}/FindQuadmath.cmake
#  COPYONLY)

export(EXPORT func-targets
  FILE ${CMAKE_CURRENT_BINARY_DIR}/FunCTargets.cmake
  NAMESPACE FunC::)

# Register package in the User Package Registry (TODO is that helpful?)
export(PACKAGE FunC)

# print every cmake variable
#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

