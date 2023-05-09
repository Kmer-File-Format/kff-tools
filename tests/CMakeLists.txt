# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set(CMAKE_BUILD_TYPE Debug)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Ofast -fno-omit-frame-pointer")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")

include_directories(${CMAKE_BINARY_DIR})
include_directories(${CMAKE_BINARY_DIR}/kff-cpp-api/)
include_directories(${CMAKE_BINARY_DIR}/src/)

set(SRCS
    main.cpp
    colinear_test.cpp
    encoding_test.cpp
    sequence_test.cpp
    compact_test.cpp
    skmers_test.cpp
)
    

# add the executable
add_executable(test ${SRCS})

# link libraries
target_link_libraries(test
    kff
    kfftl
    kffutils
)

target_include_directories(test
    PUBLIC
    "${PROJECT_SOURCE_DIR}/tests"
)