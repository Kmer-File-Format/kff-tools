
add_library(kffutils
    encoding.cpp
    sequences.cpp
    skmers.cpp
)

target_link_libraries(kffutils
    kff
)

target_include_directories(kffutils
    PUBLIC
    .
)
