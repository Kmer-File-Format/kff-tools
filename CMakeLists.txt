cmake_minimum_required(VERSION 3.10)

# set the project name
project(kff-tools)
project(kff-tools VERSION 0.1)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE Debug)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Ofast")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")

include_directories(${CMAKE_CURRENT_BINARY_DIR})

# declare libraries
add_subdirectory(kff-cpp-api/)

set(SRCS
    sequences.cpp
    encoding.cpp
    kfftools.cpp
    split.cpp
    merge.cpp
    translate.cpp
    outstr.cpp
    datarm.cpp
    disjoin.cpp
    validate.cpp
    instr.cpp
    bucket.cpp
    )
    
set(HEADERS
    sequences.hpp
    encoding.hpp
    CLI11.hpp
    kfftools.hpp
    split.hpp
    merge.hpp
    translate.hpp
    outstr.hpp
    datarm.hpp
    disjoin.hpp
    validate.hpp
    instr.hpp
    bucket.hpp
    )
# add the executable
add_executable(kff-tools ${SRCS} ${HEADERS})

# link libraries
target_link_libraries(kff-tools kff)
