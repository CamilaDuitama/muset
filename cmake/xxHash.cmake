# xxHash library

ExternalProject_Add(XXHASH
  PREFIX ${external_bindir}/XXHASH
  SOURCE_DIR ${external_dir}/xxHash
  SOURCE_SUBDIR "cmake_unofficial"
  CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DCMAKE_VERBOSE_MAKEFILE=ON
             -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
             -DBUILD_SHARED_LIBS=OFF
  INSTALL_COMMAND ""
  LOG_CONFIGURE ON
  LOG_BUILD ON
)

list(APPEND includes ${external_dir}/xxHash)
list(APPEND deps_libs ${external_bindir}/XXHASH/src/XXHASH-build/libxxhash.a)
list(APPEND deps XXHASH)
