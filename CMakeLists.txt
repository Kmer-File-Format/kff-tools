cmake_minimum_required(VERSION 3.19)

# set the project name
project(kff-tools)
project(kff-tools VERSION 0.1)

set(CMAKE_BUILD_TYPE Debug)

# set binary output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# declare libraries
add_subdirectory(kff-cpp-api/)
add_subdirectory(src/)
add_subdirectory(apps/)
add_subdirectory(tests/)

