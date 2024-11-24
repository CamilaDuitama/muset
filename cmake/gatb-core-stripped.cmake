# gatb-core-stripped library

string(REPLACE "," " " KMER_LIST_GATB ${KMER_LIST})

ExternalProject_Add(GATB
  PREFIX ${external_bindir}/GATB
  SOURCE_DIR ${external_dir}/gatb-core-stripped
  CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DCMAKE_VERBOSE_MAKEFILE=ON
             -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
             -DKSIZE_LIST=${KMER_LIST_GATB}
  INSTALL_COMMAND ""
  LOG_CONFIGURE ON
  LOG_BUILD ON
)

list(APPEND includes ${external_dir}/gatb-core-stripped/src
    ${external_dir}/gatb-core-stripped/thirdparty
    ${external_bindir}/GATB/src/GATB-build/include)

list(APPEND deps_libs ${external_bindir}/GATB/src/GATB-build/lib/Release/libgatbcore.a)
list(APPEND deps GATB)