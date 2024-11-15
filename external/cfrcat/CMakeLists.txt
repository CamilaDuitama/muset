cmake_minimum_required(VERSION 3.10.0)
project(cfrcat VERSION 0.0.1)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
set (CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set (CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}")

add_executable(cfrcat main.cpp)
target_include_directories(cfrcat PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_compile_definitions(cfrcat PRIVATE -DCFRCAT_VERSION="${CMAKE_PROJECT_VERSION}")