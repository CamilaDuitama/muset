#include <cstdlib>
#include <filesystem>
#include <memory>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <kmtricks/cmd.hpp>
#include <kmtricks/loop_executor.hpp>

#include <kmat_tools/cmd.h>

#include "muset_cli.h"


void kmtricks_pipeline(muset::muset_options_t muset_opt) {
    
    // set kmtricks pipeline options
    auto kmtricks_opt = std::make_shared<km::all_options>();
    kmtricks_opt->fof = muset_opt->fof; // --file
    kmtricks_opt->dir = muset_opt->out_dir; // --run-dir
    kmtricks_opt->kmer_size = muset_opt->kmer_size; // --kmer-size
    kmtricks_opt->c_ab_min = muset_opt->min_abundance; // --hard-min
    kmtricks_opt->lz4 = muset_opt->lz4; // --cpr
    kmtricks_opt->r_min = muset_opt->min_nb_absent; // --recurrence-min
    kmtricks_opt->logan = muset_opt->logan; // --logan
    kmtricks_opt->nb_threads = muset_opt->nb_threads; // --treads
    
    // other kmtricks options that might need to be set manually
    kmtricks_opt->count_format = km::COUNT_FORMAT::KMER; // --mode
    kmtricks_opt->mode = km::MODE::COUNT;
    kmtricks_opt->format = km::FORMAT::BIN;
    kmtricks_opt->m_ab_min = 1;
    kmtricks_opt->until = km::COMMAND::ALL;
    kmtricks_opt->minim_size = 10; // --minimizer-size
    kmtricks_opt->restrict_to = 1.0; // --restrict-to
    kmtricks_opt->focus = 0.5; // --focus
    kmtricks_opt->bloom_size = 10000000;
    kmtricks_opt->bwidth = 2;
    kmtricks_opt->out_format = km::OUT_FORMAT::HOWDE;
    kmtricks_opt->verbosity = "info";

    km::const_loop_executor<0, KMER_N>::exec<km::main_all>(kmtricks_opt->kmer_size, kmtricks_opt);
}

void kmat_filter(muset::muset_options_t muset_opt) {

    auto filter_opt = std::make_shared<kmat::filter_options>();

    filter_opt->output = muset_opt->filtered_matrix;
    
    filter_opt->min_abundance = muset_opt->min_abundance;
    
    filter_opt->min_frac_absent = muset_opt->min_frac_absent;
    filter_opt->min_frac_present = muset_opt->min_frac_present;
    filter_opt->min_nb_absent_set = muset_opt->min_nb_absent_set;
    filter_opt->min_nb_absent = muset_opt->min_nb_absent;
    filter_opt->min_nb_present_set = muset_opt->min_nb_present_set;
    filter_opt->min_nb_present = muset_opt->min_nb_present;

    filter_opt->keep_tmp = muset_opt->keep_tmp;
    filter_opt->lz4 = muset_opt->lz4;
    filter_opt->kmer_size = muset_opt->kmer_size; // not actually needed
    
    filter_opt->nb_threads = muset_opt->nb_threads;

    (filter_opt->inputs).push_back(muset_opt->out_dir);

    kmat::main_filter(filter_opt);
}


int main(int argc, char* argv[])
{
    muset::musetCli cli("muset", "a pipeline for building an abundance unitig matrix from a list of FASTA/FASTQ files.", PROJECT_VER, "");

    try
    {
        auto muset_opt = cli.parse(argc, argv);
        muset_opt->sanity_check();
        
        if (!muset_opt->min_utg_len_set) {
            muset_opt->min_utg_len = 2 * muset_opt->kmer_size - 1;
        }

        // set_verbosity_level(muset_opt->verbosity);
        auto cerr_logger = spdlog::stderr_color_mt("kmtricks");
        cerr_logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");
        spdlog::set_default_logger(cerr_logger);
        
        // auto kmer_size = muset_opt->kmer_size;
        // std::cout << muset_opt->display() << std::endl;

        if(!muset_opt->fof.empty()) {
            
            kmtricks_pipeline(muset_opt);

            muset_opt->filtered_matrix = muset_opt->out_dir/"filtered_matrix.mat";
            kmat_filter(muset_opt);

        }
    }
    catch (const km::km_exception& e) {
        spdlog::error("{} - {}", e.get_name(), e.get_msg());
        std::exit(EXIT_FAILURE);
    }
    catch (const std::exception& e) {
        spdlog::error(e.what());
        std::exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}
