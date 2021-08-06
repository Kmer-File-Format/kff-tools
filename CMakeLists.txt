cmake_minimum_required(VERSION 3.10)

# set the project name
project(kff-tools)
project(kff-tools VERSION 0.1)

set(CMAKE_BUILD_TYPE Debug)

# declare libraries
add_subdirectory(kff-cpp-api/)
add_subdirectory(src/)
add_subdirectory(tests/)

