# lz4 library

ExternalProject_Add(lz4
    PREFIX ${external_bindir}/lz4
    SOURCE_DIR ${external_dir}/lz4
    SOURCE_SUBDIR build/cmake
    CMAKE_ARGS -DBUILD_STATIC_LIBS=ON -DBUILD_SHARED_LIBS=OFF
               -DCMAKE_BUILD_TYPE=Release -DCMAKE_VERBOSE_MAKEFILE=ON
               -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    INSTALL_COMMAND ""
    LOG_CONFIGURE ON
    LOG_BUILD ON
)

list(APPEND includes ${external_dir}/lz4/lib)
list(APPEND deps_libs ${external_bindir}/lz4/src/lz4-build/liblz4.a)
list(APPEND deps lz4)
