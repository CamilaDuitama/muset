#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include "cli.h"
#include "kmtricks.h"

int main(int argc, char* argv[])
{

  muset::musetCli cli(
    "muset",
    "a pipeline for building an abundance unitig matrix from a list of FASTA/FASTQ files.",
    "v0.4.1",
    ""
  );

  try
  {
    auto options = cli.parse(argc, argv);
    options->sanity_check();

    if (!options->min_utg_len_set) {
        options->min_utg_len = 2 * options->kmer_size - 1;
    }

    // set_verbosity_level(options->verbosity);
    auto cerr_logger = spdlog::stderr_color_mt("kmtricks");
    cerr_logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");
    spdlog::set_default_logger(cerr_logger);
    
    auto kmer_size = options->kmer_size;
    // km::const_loop_executor<0, KMER_N>::exec<main_all>(kmer_size, options);
    std::cout << options->display() << std::endl;
  }
  catch (const km::km_exception& e)
  {
    spdlog::error("{} - {}", e.get_name(), e.get_msg());
    exit(EXIT_FAILURE);
  }
  catch (const std::exception& e)
  {
    spdlog::error(e.what());
    exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}