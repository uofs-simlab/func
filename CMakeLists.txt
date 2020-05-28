cmake_minimum_required (VERSION 3.1)
# SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build,
# options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
# RelWithDebInfo MinSizeRel.")

project(FunC)

set(FunC_VERSION_MAJOR 0)
set(FunC_START_VERSION_MINOR 1)

set (CMAKE_CXX_STANDARD 11)

LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

find_package(Boost COMPONENTS
        # system
        REQUIRED)

message(STATUS "Found boost incl: ${Boost_INCLUDE_DIR}")
message(STATUS "Boost libs to link: ${Boost_LIBRARIES}")
message(STATUS "Found boost lib: ${Boost_LIBRARY_DIRS}")

# Find Armadillo for the degree 4-7 polynomial interpolation and Pade tables
option(USE_ARMADILLO "Use Armadillo library" ON)
find_package(Armadillo COMPONENTS)

if(USE_ARMADILLO)
  message(STATUS "Found Armadillo incl: ${ARMADILLO_INCLUDE_DIR}")
  message(STATUS "Armadillo libs to link: ${ARMADILLO_LIBRARIES}")
  message(STATUS "Found Armadillo lib: ${ARMADILLO_LIBRARY_DIRS}")
endif()

# Search for quadmath (fatal error if requested but not found)
option(USE_QUADMATH "Use quadmath library" ON)
find_package(Quadmath)

# Add an option for turning on argument recording
option(FUNC_RECORD "Enable recording of arguments for DirectEvaluation" OFF)
message(STATUS "Argument recording is set to ${FUNC_RECORD} for DirectEvaluation")

# Add an option for removing recording arguments outside table ranges for failure-proof tables
option(NDEBUG "Disable recording arguments outside of failure-proof table ranges" OFF)
message(STATUS "Recording arguments outside of failure-proof table ranges is set to ${NDEBUG}")

# Primary source subdirectory
add_subdirectory("src")

# Examples subdirectory
add_subdirectory("examples")
