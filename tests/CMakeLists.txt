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
#add_executable(test main.cpp sequence_test.cpp)

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