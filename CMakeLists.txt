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


# Search for quadmath (fatal error if requested but not found)
option(USE_QUADMATH "Use quadmath library" ON)
find_package(Quadmath)


# Primary source subdirectory
add_subdirectory("src")

# Examples subdirectory
add_subdirectory("examples")
