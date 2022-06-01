cmake_minimum_required(VERSION 3.1)

# Default to "Release" build
SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
# options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) DEBUG RELEASE
# RELWITHDEBINFO MINSIZEREL DEV.")

project(
  FunC
  VERSION 2.0
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)

LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

option(USE_BOOST "Use Boost library" ON)

find_package(Boost COMPONENTS QUIET)

if(NOT USE_BOOST)
  message(STATUS "Not using Boost")
endif()

if(USE_BOOST AND Boost_FOUND)
  message(STATUS "Version of Boost found: ${Boost_VERSION}")
  message(STATUS "Found boost incl: ${Boost_INCLUDE_DIR}")
  message(STATUS "Boost libs to link: ${Boost_LIBRARIES}")
  message(STATUS "Found boost lib: ${Boost_LIBRARY_DIRS}")
endif()

if(USE_BOOST AND (NOT Boost_FOUND OR Boost_VERSION VERSION_LESS 1.71.0))
  set(USE_BOOST OFF)
  message(WARNING "CMake either could not find Boost or the copy of Boost found has version less than 1.71.0. You may still read table data from an existing json file but table generation is severly limited. If you wish to generate any lookup tables then you need a newer version of Boost. Silence this warning by setting -DUSE_BOOST=OFF")
endif()

# Find Armadillo for the degree 4-7 polynomial interpolation and Pade tables
option(USE_ARMADILLO "Use Armadillo library" ON)

find_package(Armadillo COMPONENTS QUIET)

if(NOT USE_ARMADILLO)
  message(STATUS "Not using Armadillo")
endif()

if(USE_ARMADILLO AND ARMADILLO_FOUND)
  message(STATUS "Found Armadillo incl: ${ARMADILLO_INCLUDE_DIR}")
  message(STATUS "Armadillo libs to link: ${ARMADILLO_LIBRARIES}")
  message(STATUS "Found Armadillo lib: ${ARMADILLO_LIBRARY_DIRS}")
endif()

if(USE_ARMADILLO AND NOT ARMADILLO_FOUND)
  set(USE_ARMADILLO OFF)
  message(WARNING "Armadillo was not found. You can read higher order tables or Pade tables from existing json files, but these table types cannot be generated. If you wish to generate higher order tables or pade tables then reconfigure with Armadillo installed. Silence this warning by setting -DUSE_ARMADILLO=OFF")
endif()

# Use quadmath if possible. Don't make a fuss if it's not there
# Could be reconfigured if the user really doesn't want quadmath
# even if they have it
find_package(Quadmath QUIET)
option(USE_QUADMATH "Use quadmath library" ON)

if(USE_QUADMATH AND NOT QUADMATH_FOUND)
  # we'll use long double instead, but tell the user what they're missing if they care
  set(USE_QUADMATH OFF)
  message(WARNING "Could not find quadmath. If you wish to use higher accuracy error estimates with the LookupTableGenerator then reconfigure to use the gnu or intel c++ compilers: ie set -DCMAKE_CXX_COMPILER=(g++-XX|icpc) (g++ on MacOS might be an alias for clang++). Silence this warning by setting -DUSE_QUADMATH=OFF")
endif()

# openmp is really not required
find_package(OpenMP QUIET)
if(OpenMP_FOUND)
  if(APPLE AND ${CMAKE_VERSION} VERSION_LESS 3.12.0)
    message(WARNING "Must have at least CMake version 3.12 to add openmp to FunC on MacOS")
    option(USE_OPENMP "Use OpenMP" OFF)
  else()
    message(STATUS "Using OpenMP for threadsafety")
    option(USE_OPENMP "Use OpenMP" ON)
  endif()
else()
  # use std::lib instead
  option(USE_OPENMP "Use OpenMP" OFF)
endif()


# Add an option for turning on argument recording
# TODO FunC is header only now, we should remove this, right?
# option(FUNC_RECORD "Enable recording of arguments for DirectEvaluation" OFF)
# message(STATUS "Argument recording is set to ${FUNC_RECORD} for DirectEvaluation")

# Primary source subdirectory
add_subdirectory("src")

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
