# kmtricks
set(KMER_LIST "32 64")
set(MAX_C 4294967295)

ExternalProject_Add(KMTRICKS
  PREFIX ${external_bindir}/KMTRICKS
  SOURCE_DIR ${external_dir}/kmtricks-logan
  CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DCMAKE_VERBOSE_MAKEFILE=ON
             -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
             -DSTATIC=OFF -DNATIVE=${ARCH_NATIVE}
             -DKMER_LIST=${KMER_LIST} -DMAX_C=${MAX_C} -DWITH_MODULES=ON
  INSTALL_COMMAND ${CMAKE_COMMAND} -E copy ${external_bindir}/KMTRICKS/src/KMTRICKS-build/src/kmtricks ${PROJECT_SOURCE_DIR}/bin/muset-kmtricks
  COMMAND make clean
  LOG_CONFIGURE ON
  LOG_BUILD ON
)

list(APPEND includes ${external_dir}/kmtricks-logan/include)
list(APPEND deps KMTRICKS)

list(APPEND includes ${external_dir}/kmtricks-logan/thirdparty/bcli/include)
list(APPEND includes ${external_dir}/kmtricks-logan/thirdparty/spdlog/include)
list(APPEND deps_libs ${external_bindir}/KMTRICKS/src/KMTRICKS-build/thirdparty/SPDLOG/src/SPDLOG-build/libspdlog.a)