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

set(SRCS
    bucket.cpp
    compact.cpp
    datarm.cpp
    disjoin.cpp
    encoding.cpp
    index.cpp
    instr.cpp
    kfftools.cpp
    merge.cpp
    outstr.cpp
    sequences.cpp
    shuffle.cpp
    sort.cpp
    split.cpp
    translate.cpp
    validate.cpp
    )
    
set(HEADERS
    CLI11.hpp
    bucket.hpp
    compact.hpp
    datarm.hpp
    disjoin.hpp
    encoding.hpp
    index.hpp
    instr.hpp
    kfftools.hpp
    merge.hpp
    outstr.hpp
    sequences.hpp
    shuffle.hpp
    sort.hpp
    split.hpp
    translate.hpp
    validate.hpp
    )

# add the executable
add_executable(kff-tools ${SRCS} ${HEADERS})

# link libraries
find_package(OpenMP)
target_link_libraries(kff-tools PUBLIC kff OpenMP::OpenMP_CXX)
