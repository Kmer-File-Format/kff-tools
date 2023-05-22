add_library(kfftl
    bucket.cpp
    compact.cpp
    datarm.cpp
    disjoin.cpp
    index.cpp
    instr.cpp
    merge.cpp
    outstr.cpp
    shuffle.cpp
    sort.cpp
    split.cpp
    translate.cpp
    validate.cpp
)
    
target_include_directories(kfftl
    PUBLIC
    ./include/
    "${PROJECT_SOURCE_DIR}/extern/CLI11/include/"
    "${PROJECT_SOURCE_DIR}/kff-cpp-api/"
)


# link libraries
find_package(OpenMP)
target_link_libraries(kfftl
    PUBLIC
    kff
    kffutils
    OpenMP::OpenMP_CXX
)
