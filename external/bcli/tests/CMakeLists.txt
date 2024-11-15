file(GLOB TEST_FILES "test_*.cpp")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/tests/)
add_executable(bcli_tests "googletest_main.cpp" ${TEST_FILES})
target_include_directories(bcli_tests PUBLIC ${PROJECT_SOURCE_DIR}/include ${GOOGLE_TEST_INCLUDE})
target_link_directories(bcli_tests PUBLIC ${GOOGLE_TEST_LIB})
target_link_libraries(bcli_tests gtest gtest_main)

add_test(
    NAME bcli_tests
    COMMAND sh -c "cd ${PROJECT_SOURCE_DIR}/tests/ ; ./bcli_tests --verbose"
)
