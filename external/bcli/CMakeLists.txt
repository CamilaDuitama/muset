cmake_minimum_required(VERSION 3.10.0)
project(bcli VERSION 0.0.1)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG")

option(COMPILE_TESTS "Compile bcli tests." OFF)
option(COMPILE_EXAMPLES "Compile bcli examples." OFF)

add_library(bcli INTERFACE)
target_include_directories(bcli INTERFACE "${PROJECT_SOURCE_DIR}/include")

if (COMPILE_TESTS)
  add_subdirectory(thirdparty/googletest)
  set(GOOGLE_TEST_INCLUDE ${PROJECT_SOURCE_DIR}/thirdparty/googletest/googletest/include)
  set(GOOGLE_TEST_LIB ${PROJECT_BINARY_DIR}/lib)
  enable_testing()
  add_subdirectory(tests)
endif()

if (COMPILE_EXAMPLES)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/ex_bin")
  add_subdirectory(examples)
endif()


