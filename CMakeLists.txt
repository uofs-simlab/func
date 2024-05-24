cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

# Default to "Release" build
set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
# options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) DEBUG RELEASE
# RELWITHDEBINFO MINSIZEREL.")

project(
  FunC
  VERSION 2.1
  DESCRIPTION "Build better lookup tables"
  LANGUAGES CXX)


include(CheckCXXCompilerFlag)
unset(FUNC_USE_CXX17 CACHE)
CHECK_CXX_COMPILER_FLAG(-std=c++17 FUNC_USE_CXX17)

if(FUNC_USE_CXX17)
set(CMAKE_CXX_STANDARD 17)
else()
set(CMAKE_CXX_STANDARD 14)
endif()

macro(print_variable var)
  message(STATUS "${var}: ${${var}}")
endmacro()

# cmake will look in our cmake directory for FindXXX.cmake files
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

# ==== Search for dependencies ====
option(FUNC_USE_BOOST "Use Boost library" ON)
set(MIN_BOOST_VERSION 1.71.0)
# TODO the following line may cause the build to fail if Boost doesn't exist on the user's computer.
# It seems like CMake's FindBoost is the problem?
find_package(Boost ${MIN_BOOST_VERSION})

if(FUNC_USE_BOOST)
  if(Boost_FOUND AND NOT (Boost_VERSION VERSION_LESS ${MIN_BOOST_VERSION}))
    message(STATUS "Found boost incl: ${Boost_INCLUDE_DIR}")
  else()
    message(FATAL_ERROR "CMake either could not find Boost or Boost was found with version less than 1.71.0. Without Boost, FunC can't generate most LUT types unless they are read from an existing JSON file. Further, FunC cannot compute a LUT's error without Boost. To use FunC without Boost, run CMake with the option: -DFUNC_USE_BOOST=OFF")
  endif()
else()
  message(STATUS "Not using Boost")
endif()


# Find Armadillo for the degree 4-7 polynomial interpolation and Pade tables
option(FUNC_USE_ARMADILLO "Use Armadillo library" ON)
find_package(Armadillo)

if(FUNC_USE_ARMADILLO)
  if(ARMADILLO_FOUND)
    message(STATUS "Armadillo includes: ${ARMADILLO_INCLUDE_DIRS}")
    message(STATUS "Armadillo libs to link: ${ARMADILLO_LIBRARIES}")
  else()
    message(FATAL_ERROR "Armadillo was not found. Without Armadillo, FunC cannot generate ChebyInterpTables or PadeTable unless they are read from an existing JSON file. To use FunC without Armadillo, run CMake with the option: -DFUNC_USE_ARMADILLO=OFF")
  endif()
else()
  message(STATUS "Not using Armadillo")
endif()

option(FUNC_USE_OPENMP "Use OpenMP" ON)
find_package(OpenMP)

if(FUNC_USE_OPENMP)
  if(OpenMP_FOUND)
    message(STATUS "OpenMP lib to incl: ${OpenMP_CXX_LIBRARIES}")
    set(OPENMP_CXX_LIBRARIES OpenMP::OpenMP_CXX)
  else()
    message(FATAL_ERROR "OpenMP was not found. Without OpenMP, FunC cannot parallelize the processes of building a LUT or checking its error. To use FunC without OpenMP, run CMake with the option: -DFUNC_USE_OPENMP=OFF")
  endif()
else()
  message(STATUS "Not using OpenMP")
endif()

# Any conversions from json type should be done explicitly
set(JSON_ImplicitConversions OFF)

# ==== Generate Source Code ====
# TYPE_PAIRS contains every type we want to use out LUTs with.
# TODO can the average user change TYPE_PAIRS without wanting to cry?
if(NOT DEFINED TYPE_PAIRS)
  message(STATUS "Variable TYPE_PAIRS was not set by user")
  set(TYPE_PAIRS "double,double;float,float;long double,long double")
endif()

# Every classname with whose source we want compiled into libfunc.so for each type_pair
# This is the only file that showed material speedup
set(CLASSNAMES
  LookupTableFactory # .hpp file must contain a line with FUNC_DECLARE_TEMPLATE_AS_EXTERN(classname)
  LookupTableComparator
)

# cleanup any source-code that may have been generated in previous cmake runs
set(func_src_dir ${CMAKE_CURRENT_BINARY_DIR}/src)
file(GLOB old_func_src ${func_src_dir}/*.cpp)
if(NOT "${old_func_src}" STREQUAL "")
  file(REMOVE ${old_func_src})
endif()

set(FUNC_DECLARE_TEMPLATE_AS_EXTERN "#define FUNC_DECLARE_TEMPLATE_AS_EXTERN(classname)\\")

foreach(PAIR IN LISTS TYPE_PAIRS)
  set(FUNC_DECLARE_TEMPLATE_AS_EXTERN "${FUNC_DECLARE_TEMPLATE_AS_EXTERN}\n  extern template class classname<${PAIR}>;\\")

  string(REPLACE "," "_to_" FORMAT_PAIR "${PAIR}")
  foreach(CLASSNAME IN LISTS CLASSNAMES)
    configure_file(${PROJECT_SOURCE_DIR}/cmake/src_template.cpp.cmake ${func_src_dir}/${CLASSNAME}_${FORMAT_PAIR}.cpp @ONLY)
  endforeach()
endforeach()

# declare the extern template types and let c++ know what dependencies we found
configure_file(${PROJECT_SOURCE_DIR}/cmake/config.hpp.cmake ${PROJECT_SOURCE_DIR}/include/config.hpp @ONLY)

# Primary source subdirectory
add_subdirectory("include")

# Examples subdirectory, if enabled
option(BUILD_EXAMPLES "Build the example codes" ON)
if(BUILD_EXAMPLES AND CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)
  add_subdirectory("examples")
endif()


# Summarize compilation flags for quick review on configuration
macro(print_variable var)
  message(STATUS "${var}: ${${var}}")
endmacro()

message(STATUS "CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message(STATUS "--- Compile and link flags ------------------")
message(STATUS "CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
string(TOLOWER "${CMAKE_BUILD_TYPE}" build_lower)
if(build_lower MATCHES release)
  print_variable(CMAKE_CXX_FLAGS_RELEASE)
elseif(build_lower MATCHES debug)
  print_variable(CMAKE_CXX_FLAGS_DEBUG)
elseif(build_lower MATCHES relwithdebinfo)
  print_variable(CMAKE_CXX_FLAGS_RELWITHDEBINFO)
elseif(build_lower MATCHES minsizerel)
  print_variable(CMAKE_CXX_FLAGS_MINSIZEREL)
endif()

# Print every variable for debugging purposes
#get_cmake_property(_variableNames VARIABLES)
#list (SORT _variableNames)
#foreach (_variableName ${_variableNames})
#    message(STATUS "${_variableName}=${${_variableName}}")
#endforeach()

message(STATUS "---------------------------------------------")

