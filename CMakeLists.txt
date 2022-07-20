cmake_minimum_required(VERSION 3.1)

# Default to "Release" build
SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
# options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) DEBUG RELEASE
# RELWITHDEBINFO MINSIZEREL DEV.")

project(
  FunC
  VERSION 2.0
  DESCRIPTION "Build better lookup tables"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

# ==== Search for dependencies ====
option(FUNC_USE_BOOST "Use Boost library" ON)
find_package(Boost 1.71.0)

if(FUNC_USE_BOOST)
  if(Boost_FOUND AND NOT (Boost_VERSION VERSION_LESS 1.71.0))
    message(STATUS "Found boost incl: ${Boost_INCLUDE_DIR}")
  else()
    set(FUNC_USE_BOOST OFF)
    message(WARNING "CMake either could not find Boost or the copy of Boost found has version less than 1.71.0. You may still read table data from an existing json file but table generation is severly limited. If you wish to generate any lookup tables then you need a newer version of Boost. Silence this warning by setting -DUSE_BOOST=OFF")
  endif()
else()
  message(STATUS "Not using Boost")
endif()


# Find Armadillo for the degree 4-7 polynomial interpolation and Pade tables
option(FUNC_USE_ARMADILLO "Use Armadillo library" ON)
find_package(Armadillo)

if(FUNC_USE_ARMADILLO)
  if(ARMADILLO_FOUND)
    message(STATUS "Found Armadillo incl: ${ARMADILLO_INCLUDE_DIR}")
    message(STATUS "Armadillo libs to link: ${ARMADILLO_LIBRARIES}")
    message(STATUS "Found Armadillo lib: ${ARMADILLO_LIBRARY_DIRS}")
  else()
    set(FUNC_USE_ARMADILLO OFF)
    message(WARNING "Armadillo was not found. You can read higher order tables or Pade tables from existing json files, but these table types cannot be generated. If you wish to generate higher order tables or pade tables then reconfigure with Armadillo installed. Silence this warning by setting -DUSE_ARMADILLO=OFF")
  endif()
else()
  message(STATUS "Not using Armadillo")
endif()

# use FindQuadmath.cmake so we can write find_package(Quadmath)
# This will mostly just test whether we're compiling with gcc/icpc
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
option(FUNC_USE_QUADMATH "Use quadmath library" ON)
find_package(Quadmath)

if(FUNC_USE_QUADMATH AND NOT QUADMATH_FOUND)
  # we'll use long double instead
  set(FUNC_USE_QUADMATH OFF)
  message(WARNING "Could not find quadmath. If you wish to use higher accuracy error estimates with the LookupTableGenerator then reconfigure to use the gnu or intel c++ compilers: ie set -DCMAKE_CXX_COMPILER=(g++-XX|icpc) (g++ on MacOS might be an alias for clang++). Silence this warning by setting -DUSE_QUADMATH=OFF")
endif()

option(FUNC_USE_OPENMP "Use OpenMP" ON)
find_package(OpenMP)
if(OpenMP_FOUND)
  if(APPLE AND ${CMAKE_VERSION} VERSION_LESS 3.12.0)
    message(WARNING "Must have at least CMake version 3.12 to add openmp to FunC on MacOS")
    set(FUNC_USE_OPENMP OFF)
  else()
    message(STATUS "Using OpenMP for threadsafety")
  endif()
else()
  # TODO use std::lib instead
  option(FUNC_USE_OPENMP OFF)
endif()


# ==== Generate Source Code ====
# Every type that we want to use our LUTs with. IN_TYPES and OUT_TYPES must be the same length
# TODO could the average user change these variables without wanting to cry? Does "option" support lists??
set(TYPE_PAIRS "double,double")

# Every classname with whose source we want compiled into libfunc.so for each type_pair
# This is the only file that showed material speedup
set(CLASSNAMES
  LookupTableFactory
)

# cleanup any source-code that may have been generated in previous cmake runs
set(func_src_dir ${CMAKE_CURRENT_BINARY_DIR}/src)
file(GLOB old_func_src ${func_src_dir}/*.cpp)
if(NOT "${old_func_src}" STREQUAL "")
  file(REMOVE ${old_func_src})
endif()

set(FUNC_DECLARE_TEMPLATE_AS_EXTERN "#define FUNC_DECLARE_TEMPLATE_AS_EXTERN(classname)\\")

# TODO do these configure_files need some fancy generator-ish path names?
foreach(PAIR IN LISTS TYPE_PAIRS)
  set(FUNC_DECLARE_TEMPLATE_AS_EXTERN "${FUNC_DECLARE_TEMPLATE_AS_EXTERN}\n  extern template class classname<${PAIR}>;\\")
  foreach(CLASSNAME IN LISTS CLASSNAMES)
    # TODO can we leave this automatically generated source in ${CMAKE_CURRENT_BINARY_DIR}?
    string(REPLACE "," "_to_" FORMAT_PAIR "${PAIR}")
    configure_file(${PROJECT_SOURCE_DIR}/include/src_template.cpp.cmake ${func_src_dir}/${CLASSNAME}_${FORMAT_PAIR}.cpp @ONLY)
  endforeach()
endforeach()

# declare the extern template types and let c++ know what dependencies we found
configure_file(${PROJECT_SOURCE_DIR}/include/config.hpp.cmake ${PROJECT_SOURCE_DIR}/include/config.hpp @ONLY)

# Primary source subdirectory
add_subdirectory("include")

# Examples subdirectory
add_subdirectory("examples")

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

message(STATUS "---------------------------------------------")
