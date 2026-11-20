set(SRCS
    main.cpp
    colinear_test.cpp
    encoding_test.cpp
    sequence_test.cpp
    compact_test.cpp
    skmers_test.cpp
    utils_test.cpp
    outstr_test.cpp
    blocks_test.cpp
    colinear_reference_test.cpp
    compact_key_test.cpp
    polish_reference_test.cpp
    pair_reference_test.cpp
    chain_columns_test.cpp
    greedy_reference_test.cpp
    minimizer_reference_test.cpp
)

# Unit tests ("test" is a target name reserved by CTest)
add_executable(kff-tools-tests ${SRCS})
set_target_properties(kff-tools-tests PROPERTIES OUTPUT_NAME tests)

# link libraries
target_link_libraries(kff-tools-tests
    kff
    kfftl
    kffutils
)

target_include_directories(kff-tools-tests
    PUBLIC
    "${PROJECT_SOURCE_DIR}/tests"
)

add_test(NAME unit COMMAND kff-tools-tests)

# Python functional tests (they call ./bin/kff-tools from the project root)
find_package(Python3 COMPONENTS Interpreter)
if (Python3_FOUND)
    add_test(NAME functional
        COMMAND ${Python3_EXECUTABLE} -m unittest discover -s tests
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    )
endif()
