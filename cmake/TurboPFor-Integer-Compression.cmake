# TurboPFor-Integer-Compression library

ExternalProject_Add(TURBOP
  PREFIX ${external_bindir}/TURBOP
  SOURCE_DIR ${external_dir}/TurboPFor-Integer-Compression
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ${CMAKE_COMMAND} -E env
  make libic.a
  COMMAND ${CMAKE_COMMAND} -E copy "<SOURCE_DIR>/libic.a" ${external_bindir}/TURBOP
  INSTALL_COMMAND make clean
  LOG_CONFIGURE ON
  LOG_BUILD ON
)

list(APPEND includes ${external_dir}/TurboPFor-Integer-Compression/include)
list(APPEND deps_libs ${external_bindir}/TURBOP/libic.a)
list(APPEND deps TURBOP)
