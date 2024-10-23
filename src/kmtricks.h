#define KMER_LIST 32, 64
#define KMER_N 2
#define DMAX_C 4294967295
#define WITH_KM_IO

#define KMTRICKS_PUBLIC
#include <kmtricks/kmer.hpp>
// #include <kmtricks/minimizer.hpp>
#include <kmtricks/repartition.hpp>
// #include <kmtricks/histogram.hpp>
// #include <kmtricks/bitmatrix.hpp>
// #include <kmtricks/hash.hpp>
// #include <kmtricks/timer.hpp>
#include <kmtricks/itask.hpp>
// #include <kmtricks/hash.hpp>
#include <kmtricks/loop_executor.hpp>
#include <kmtricks/task_pool.hpp>

// require lz4
#include <kmtricks/io/matrix_file.hpp>