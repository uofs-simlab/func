cmake_minimum_required (VERSION 3.1)

# Default to "Release" build
SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
# options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) DEBUG RELEASE
# RELWITHDEBINFO MINSIZEREL.")

project(FunC)

set(FunC_VERSION_MAJOR 0)
set(FunC_START_VERSION_MINOR 1)

set (CMAKE_CXX_STANDARD 11)

LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

find_package(Boost COMPONENTS
        # system
        REQUIRED)

message(STATUS "Version of Boost found: ${Boost_VERSION}")
message(STATUS "Found boost incl: ${Boost_INCLUDE_DIR}")
message(STATUS "Boost libs to link: ${Boost_LIBRARIES}")
message(STATUS "Found boost lib: ${Boost_LIBRARY_DIRS}")

if(Boost_VERSION VERSION_LESS 1.71.0)
  set(USE_BOOST_AUTODIFF OFF)
  message(WARNING "The version of Boost found is less than 1.71.0
If you wish to use Taylor, Pade, and Hermite tables then you will need
a newer version of Boost.")
else()
  set(USE_BOOST_AUTODIFF ON)
endif()

# Find Armadillo for the degree 4-7 polynomial interpolation and Pade tables
option(USE_ARMADILLO "Use Armadillo library" OFF)

find_package(Armadillo COMPONENTS)

if(ARMADILLO_FOUND)
  if(NOT USE_ARMADILLO)
    message(WARNING "Armadillo found but not requested.
If you wish to use Armadillo (NEEDED FOR HIGHER ORDER AND PADE TABLES),
 reconfigure with -DUSE_ARMADILLO=ON")
  endif()
  message(STATUS "Found Armadillo incl: ${ARMADILLO_INCLUDE_DIR}")
  message(STATUS "Armadillo libs to link: ${ARMADILLO_LIBRARIES}")
  message(STATUS "Found Armadillo lib: ${ARMADILLO_LIBRARY_DIRS}")
endif()

# Use quadmath if possible. Don't make a fuss if it's not there
# Could be reconfigured if the user really doesn't want quadmath
# even if they have it
find_package(Quadmath)
if(QUADMATH_FOUND)  
  option(USE_QUADMATH "Use quadmath library" ON)
else()
  # we'll just use long double instead
  option(USE_QUADMATH "Use quadmath library" OFF)
endif()

# Add an option for turning on argument recording
option(FUNC_RECORD "Enable recording of arguments for DirectEvaluation" OFF)
message(STATUS "Argument recording is set to ${FUNC_RECORD} for DirectEvaluation")

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
