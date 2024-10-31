cmake_minimum_required(VERSION 3.15)
project(kmat_tools LANGUAGES C CXX)

############################################################
## Required C++ standard

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(cmake_path "${PROJECT_SOURCE_DIR}/cmake")

############################################################
## Options

option(ARCH_NATIVE "Add -march=native compiler flag. Disabled for conda builds." ON)
option(CONDA_BUILD "Build inside conda env." OFF)

############################################################
## Prevent in-source build

if(${PROJECT_SOURCE_DIR} STREQUAL ${PROJECT_BINARY_DIR})
    message(FATAL_ERROR "In-source build not allowed. Create a \"build\" directory and run CMake from there.")
endif()

############################################################
## Default build type set to "Release"

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()
message(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})

############################################################
## Disable -march=native flag for conda builds

if (CONDA_BUILD)
  set(ARCH_NATIVE OFF)
endif()

############################################################
## Display details on host system

message(STATUS "Build Info:")
message(STATUS "  Host System: ${CMAKE_HOST_SYSTEM}")
message(STATUS "  Host System Name: ${CMAKE_HOST_SYSTEM_NAME}")
message(STATUS "  Host System Processor: ${CMAKE_HOST_SYSTEM_PROCESSOR}")
message(STATUS "  Host System Version: ${CMAKE_HOST_SYSTEM_VERSION}")

############################################################
## Compiler flags

set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-ignored-qualifiers -Wno-char-subscripts -Wno-maybe-uninitialized -Wno-unused-function -pedantic")

if (UNIX)
  if (CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "amd64|x86_64")
    # flags for PTHash
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2 -msse4.2")
  endif()
endif()

if (ARCH_NATIVE)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
endif()

message(STATUS "  CXX Compiler flags: ${CMAKE_CXX_FLAGS}")

############################################################
## Dependencies

include(ExternalProject)

set(external_dir ${PROJECT_SOURCE_DIR}/external)
set(external_bindir ${PROJECT_BINARY_DIR}/external)

set(includes "")
set(deps_libs "")
set(deps "")

# kmtricks
include(${cmake_path}/kmtricks.cmake)

# fmt
include(${cmake_path}/fmt.cmake)

# lz4
include(${cmake_path}/lz4.cmake)

# sshash
include(${cmake_path}/sshash.cmake)

# zlib
find_package(ZLIB REQUIRED)
list(APPEND deps_libs ZLIB::ZLIB)

# threads
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)
list(APPEND deps_libs Threads::Threads)

############################################################
## Set executable output directory

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
message(STATUS "Executables output directory: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

############################################################
# kmat_tools executable

set(kmat_tools_sources
  src/km_basic_filter.cpp
  src/km_diff.cpp
  src/km_fafmt.cpp
  src/km_fasta.cpp
  src/km_ktfilter.cpp
  src/km_merge.cpp
  src/km_reverse.cpp
  src/km_select.cpp
  src/km_tools.cpp
  src/km_unitig.cpp
  src/km_convert.cpp
)

add_executable(kmat_tools ${kmat_tools_sources})

target_include_directories(kmat_tools PRIVATE ${includes})
target_link_libraries(kmat_tools ${deps_libs})
add_dependencies(kmat_tools ${deps})

############################################################
# muset cpp executable

set(muset_sources
  src/muset.cpp
  src/cli.cpp
)

add_executable(muset-cpp ${muset_sources})

target_include_directories(muset-cpp PRIVATE ${includes})
target_link_libraries(muset-cpp ${deps_libs})
add_dependencies(muset-cpp ${deps})
