# spdlog library

ExternalProject_Add(SPDLOG
  PREFIX ${external_bindir}/SPDLOG
  SOURCE_DIR ${external_dir}/spdlog
  CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DCMAKE_VERBOSE_MAKEFILE=ON
             -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
             -DSPDLOG_BUILD_SHARED=OFF
  INSTALL_COMMAND ""
  LOG_CONFIGURE ON
  LOG_BUILD ON
)

list(APPEND includes ${external_dir}/spdlog/include)
list(APPEND deps SPDLOG)
