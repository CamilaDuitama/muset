# kff-cpp-api library

ExternalProject_Add(KFF
  PREFIX ${external_bindir}/KFF
  SOURCE_DIR ${external_dir}/kff-cpp-api
  CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DCMAKE_VERBOSE_MAKEFILE=ON
             -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
  INSTALL_COMMAND ""
  LOG_CONFIGURE ON
  LOG_BUILD ON
)

list(APPEND includes ${external_bindir}/KFF/src/KFF-build)
list(APPEND deps_libs ${external_bindir}/KFF/src/KFF-build/libkff.a)
list(APPEND deps KFF)
