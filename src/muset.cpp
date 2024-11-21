#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <sstream>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <kmtricks/cmd.hpp>
#include <kmtricks/loop_executor.hpp>

#include <kmat_tools/cmd/fafmt.h>
#include <kmat_tools/cmd/fasta.h>
#include <kmat_tools/cmd/filter.h>
#include <kmat_tools/cmd/unitig.h>
#include <kmat_tools/matrix.h>

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

void kmat_filter(muset::muset_options_t muset_opt, bool kmtricks_input) {

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

    (filter_opt->inputs).push_back(kmtricks_input ? muset_opt->out_dir : muset_opt->in_matrix);

    kmat::main_filter(filter_opt);
}

void kmat_fasta(muset::muset_options_t muset_opt) {

    auto fasta_opt = std::make_shared<kmat::fasta_options>();
    fasta_opt->output = muset_opt->filtered_kmers;
    (fasta_opt->inputs).push_back(muset_opt->filtered_matrix);

    kmat::main_fasta(fasta_opt);
}

void ggcat(muset::muset_options_t muset_opt) {

    auto fasta_opt = std::make_shared<kmat::fasta_options>();
    fasta_opt->output = muset_opt->filtered_kmers;
    (fasta_opt->inputs).push_back(muset_opt->filtered_matrix);

    std::string ggcat_logfile = muset_opt->out_dir/"ggcat.log";
    std::string ggcat_cmd = fmt::format("ggcat build -j {} -k {} -s 1 -o {} {} &> {}",
        muset_opt->nb_threads, muset_opt->kmer_size,
        (muset_opt->unitigs).c_str(),
        (muset_opt->filtered_kmers).c_str(),
        ggcat_logfile);

    spdlog::debug(fmt::format("Running command: {}", ggcat_cmd));    
    auto ret = std::system(ggcat_cmd.c_str());
    if(ret != 0) {
        throw std::runtime_error(fmt::format("Command failed: {}\nSee log at {}", ggcat_cmd));
    }
}

void kmat_fafmt(muset::muset_options_t muset_opt) {

    auto fafmt_opt = std::make_shared<kmat::fafmt_options>();
    fafmt_opt->output = muset_opt->filtered_unitigs;
    fafmt_opt->min_length = muset_opt->min_utg_len;
    (fafmt_opt->inputs).push_back(muset_opt->unitigs);

    kmat::main_fafmt(fafmt_opt);
}

void kmat_unitig(muset::muset_options_t muset_opt) {

    auto unitig_opt = std::make_shared<kmat::unitig_options>();

    unitig_opt->kmer_size = muset_opt->kmer_size;
    unitig_opt->mini_size = muset_opt->mini_size;
    unitig_opt->prefix = muset_opt->unitig_prefix;
    unitig_opt->min_frac = muset_opt->min_utg_frac;
    unitig_opt->write_seq = muset_opt->write_utg_seq;
    unitig_opt->write_frac_matrix = muset_opt->write_frac_matrix;
    unitig_opt->nb_threads = muset_opt->nb_threads;

    (unitig_opt->inputs).push_back(muset_opt->filtered_unitigs);
    (unitig_opt->inputs).push_back(muset_opt->filtered_matrix);

    kmat::main_unitig(unitig_opt);
}

int main(int argc, char* argv[])
{
    muset::musetCli cli("muset", "a pipeline for building an abundance unitig matrix from a list of FASTA/FASTQ files.", PROJECT_VER, "");

    auto cerr_logger = spdlog::stderr_color_mt("muset");
    cerr_logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");
    spdlog::set_default_logger(cerr_logger);

    try
    {
        // check ggcat dependency
        if (std::system("ggcat --version &>/dev/null") != 0) {
            throw std::runtime_error("ggcat: command not found");
        }

        auto muset_opt = cli.parse(argc, argv);
        muset_opt->sanity_check();
        
        if (!muset_opt->min_utg_len_set) {
            muset_opt->min_utg_len = 2 * muset_opt->kmer_size - 1;
        }

        if(!muset_opt->fof.empty()) {
            spdlog::info("Building k-mer matrix with kmtricks");
            kmtricks_pipeline(muset_opt);
            muset_opt->kmer_matrix = muset_opt->out_dir;
        } else {
            spdlog::info(fmt::format("Using input k-mer matrix: {}", (muset_opt->in_matrix).c_str()));
            
            std::string kmer; kmat::TextMatrixReader mat(muset_opt->in_matrix);
            if (!mat.read_kmer(kmer)) {
                throw std::runtime_error("Empty input matrix");
            }
            muset_opt->kmer_size = kmer.size();
            spdlog::debug(fmt::format("input matrix k-mer size: {}", muset_opt->kmer_size));
            fs::create_directories(muset_opt->out_dir);
            muset_opt->kmer_matrix = muset_opt->in_matrix;
        }

        spdlog::info(fmt::format("Filtering k-mer matrix"));
        muset_opt->filtered_matrix = muset_opt->out_dir/"matrix.filtered.mat";
        kmat_filter(muset_opt, !muset_opt->fof.empty());

        spdlog::info(fmt::format("Writing k-mers in FASTA format"));
        muset_opt->filtered_kmers = muset_opt->out_dir/"matrix.filtered.fasta";
        kmat_fasta(muset_opt);

        spdlog::info(fmt::format("Building unitigs"));
        muset_opt->unitigs = muset_opt->out_dir/"unitigs.fa";
        ggcat(muset_opt);

        spdlog::info(fmt::format("Filtering unitigs"));
        muset_opt->filtered_unitigs = muset_opt->out_dir/"unitigs.filtered.fa";
        kmat_fafmt(muset_opt);

        spdlog::info(fmt::format("Building unitig matrix"));
        muset_opt->unitig_prefix = muset_opt->out_dir/"unitigs";
        kmat_unitig(muset_opt);
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
