cmake_minimum_required(VERSION 3.10)

# set the project name
project(kff-tools)
project(kff-tools VERSION 0.1)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set(CMAKE_BUILD_TYPE Debug)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Ofast")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")

# declare libraries
add_subdirectory(kmer_file_format/C++/)

set(SRCS
    kfftools.cpp
    split.cpp
    merge.cpp
    )
    
set(HEADERS
    CLI11.hpp
    kfftools.hpp
    split.hpp
    merge.hpp
    )
# add the executable
add_executable(kff-tools ${SRCS} ${HEADERS})

# link libraries
target_link_libraries(kff-tools kff)
