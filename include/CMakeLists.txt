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
table_types/func_tables.hpp;utils/func_utils.hpp;\
table_types/TransferFunction.hpp;table_types/Polynomial.hpp;table_types/MetaTable.hpp;\
table_types/ChebyInterpTable.hpp;\
table_types/CubicHermiteTable.hpp;\
table_types/EqSpaceInterpTable.hpp;\
table_types/LinearRawInterpTable.hpp;\
table_types/PadeTable.hpp;\
table_types/TaylorTable.hpp;\
utils/DirectEvaluation.hpp;utils/ArgumentRecord.hpp;\
LookupTableFactory.hpp;utils/LookupTableGenerator.hpp;\
utils/FailureProofTable.hpp;utils/CompositeLookupTable.hpp;\
utils/RngInterface.hpp;utils/StdRng.hpp;\
utils/LookupTableComparator.hpp;\
utils/Timer.hpp;\
")

target_include_directories(func ${FUNC_PUBLIC}
  $<INSTALL_INTERFACE:.> $<INSTALL_INTERFACE:table_types> $<INSTALL_INTERFACE:utils>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/table_types> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/utils>
  )


#target_link_libraries(func ${FUNC_PUBLIC} ${ARMADILLO_LIBRARIES} ${QUADMATH_LIBRARIES} Boost::boost)
target_link_libraries(func ${FUNC_PUBLIC} ${ARMADILLO_LIBRARIES} ${QUADMATH_LIBRARIES} Boost::boost OpenMP::OpenMP_CXX)
#target_link_libraries(func ${FUNC_PUBLIC} ${ARMADILLO_LIBRARIES} ${QUADMATH_LIBRARIES} ${Boost_LIBRARIES} ${OpenMP_CXX_LIBRARIES})

# TODO add more required compile features (it's nice to have cmake determine the minimal standard on its own)
target_compile_features(func ${FUNC_PRIVATE} cxx_range_for)

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
