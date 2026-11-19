cmake_minimum_required(VERSION 3.16...3.31)

# set the project name
project(kff-tools VERSION 0.1 LANGUAGES CXX)

# Optimised build by default
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

# Same language level and options for every target of the project
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# set binary output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

# declare libraries
add_subdirectory(kff-cpp-api/)

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall)
endif()

add_subdirectory(src/)
add_subdirectory(apps/)

enable_testing()
add_subdirectory(tests/)
