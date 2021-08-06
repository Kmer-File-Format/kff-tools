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

include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${MyProject_SOURCE_DIR}/kff-cpp-api/)
include_directories(${MyProject_SOURCE_DIR}/src/)

set(SRCS
    main.cpp
    encoding_test.cpp
    ../src/sequences.cpp
    ../src/encoding.cpp
    )
    
set(HEADERS
    ../src/sequences.hpp
    ../src/encoding.hpp
    )
# add the executable
add_executable(test ${SRCS} ${HEADERS})

# link libraries
target_link_libraries(test kff)