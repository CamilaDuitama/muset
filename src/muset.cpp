#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <sstream>

#include <fmt/format.h>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>

#include <kmtricks/cmd.hpp>
#include <kmtricks/kmdir.hpp>
#include <kmtricks/loop_executor.hpp>

#include <kmat_tools/cmd/fafmt.h>
#include <kmat_tools/cmd/fasta.h>
#include <kmat_tools/cmd/filter.h>
#include <kmat_tools/cmd/unitig.h>
#include <kmat_tools/matrix.h>
#include <kmat_tools/utils.h>

#include "muset_cli.h"


void init_logger(const fs::path &dir) {
    auto cerr_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>(); // spdlog::stderr_color_mt("muset");
    cerr_sink->set_pattern("[%Y-%m-%d %H:%M:%S] [%^%l%$] %v");

    auto now = std::chrono::system_clock::now();
    auto local_time = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss; ss << std::put_time(std::localtime(&local_time), "%Y%m%d_%H%M%S");
    auto muset_log = fmt::format("{}/muset_{}.log", dir.c_str(), ss.str());
    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(muset_log, true);
    file_sink->set_pattern("[%Y-%m-%d %H:%M:%S] [%^%l%$] %v");
    // file_sink->set_level(spdlog::level::info);

    auto combined_logger = std::make_shared<spdlog::logger>("muset", spdlog::sinks_init_list({cerr_sink,file_sink}));
    spdlog::set_default_logger(combined_logger);
}

void print_options(muset::muset_options_t opt) {

    if(!(opt->fof).empty()) { spdlog::info(fmt::format("input file (--file): {}", (opt->fof).c_str())); }
    if(!(opt->in_matrix).empty()) { spdlog::info(fmt::format("input matrix (-i): {}", (opt->in_matrix).c_str())); }
    spdlog::info(fmt::format("output directory (-o): {}", (opt->out_dir).c_str()));
    
    spdlog::info(fmt::format("k-mer size (-k): {}", opt->kmer_size));
    spdlog::info(fmt::format("minimum abundance (-a): {}", opt->min_abundance));
    spdlog::info(fmt::format("minimum unitig length (-l): {}", opt->min_utg_len));
    spdlog::info(fmt::format("minimum unitig fraction (-r): {}", opt->min_utg_frac));
    spdlog::info(fmt::format("write unitig sequence (-s): {}", opt->write_utg_seq));
    spdlog::info(fmt::format("output unitig fraction matrix (--out-frac): {}", opt->write_frac_matrix));
    spdlog::info(fmt::format("input consists of logan unitigs (--logan): {}", opt->logan));
    spdlog::info(fmt::format("minimizer size (-m): {}", opt->mini_size));

    if(opt->min_nb_absent_set) {
        spdlog::info(fmt::format("number of absent samples (-n): {}", opt->min_nb_absent));
    } else {
        spdlog::info(fmt::format("fraction of absent samples (-f): {}", opt->min_frac_absent));
    }
    if(opt->min_nb_present_set) {
        spdlog::info(fmt::format("number of present samples (-N): {}", opt->min_nb_present));
    } else {
        spdlog::info(fmt::format("fraction of present samples (-F): {}", opt->min_frac_present));
    }

    spdlog::info(fmt::format("keep temporary files (--keep-temp): {}", opt->keep_tmp));
    spdlog::info(fmt::format("threads (-t): {}", opt->nb_threads));
}

void kmtricks_pipeline(muset::muset_options_t muset_opt) {
    
    // set kmtricks pipeline options
    auto kmtricks_opt = std::make_shared<km::all_options>();
    kmtricks_opt->fof = muset_opt->fof; // --file
    kmtricks_opt->dir = muset_opt->kmer_matrix; // --run-dir
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

    (filter_opt->inputs).push_back(muset_opt->kmer_matrix);

    kmat::main_filter(filter_opt);
}

void kmat_fasta(muset::muset_options_t muset_opt) {

    auto fasta_opt = std::make_shared<kmat::fasta_options>();
    fasta_opt->output = muset_opt->filtered_kmers;
    (fasta_opt->inputs).push_back(muset_opt->filtered_matrix);

    kmat::main_fasta(fasta_opt);
}

void ggcat(muset::muset_options_t muset_opt) {
    std::string ggcat_keeptmp = muset_opt->keep_tmp ? "--keep-temp-files" : "";
    std::string ggcat_tempdir = muset_opt->out_dir/"ggcat_build_temp";
    std::string ggcat_logfile = muset_opt->out_dir/"ggcat.log";
    std::string ggcat_cmd = fmt::format("ggcat build -j {} -k {} -s 1 -o {} {} --temp-dir {} {} &> {}",
        muset_opt->nb_threads, muset_opt->kmer_size,
        (muset_opt->unitigs).c_str(),
        ggcat_keeptmp, // --keep-temp-files ?
        ggcat_tempdir, // --temp-dir
        (muset_opt->filtered_kmers).c_str(),
        ggcat_logfile);

    spdlog::debug(fmt::format("Running command: {}", ggcat_cmd));    
    auto ret = std::system(ggcat_cmd.c_str());
    if(ret != 0) {
        throw std::runtime_error(fmt::format("Command failed: {}\nSee log at {}", ggcat_cmd, ggcat_logfile));
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
    auto muset_opt = cli.parse(argc, argv);

    try
    {
        // check ggcat dependency
        if (std::system("ggcat --version &>/dev/null") != 0) {
            throw std::runtime_error("ggcat: command not found");
        }

        // check parameters consistency
        muset_opt->sanity_check();
        
        if (!muset_opt->min_utg_len_set) {
            muset_opt->min_utg_len = 2 * muset_opt->kmer_size - 1;
        }

        // create main output directory
        fs::create_directories(muset_opt->out_dir);

        // initialize the logger just before running the pipeline
        init_logger(muset_opt->out_dir);

        // print values of main options
        spdlog::info(fmt::format("Running muset {}", PROJECT_VER));
        spdlog::info("---------------");
        print_options(muset_opt);
        spdlog::info("---------------");

        // muset pipeline

        if(!muset_opt->fof.empty()) {
            // create kmtricks matrix
            spdlog::info("Building k-mer matrix with kmtricks");
            muset_opt->kmer_matrix = muset_opt->out_dir/"kmer_matrix";
            if(fs::is_directory(muset_opt->kmer_matrix)) {
                throw std::runtime_error(fmt::format("kmtricks output directory \"{}\" already exists.", (muset_opt->kmer_matrix).c_str()));
            }
            kmtricks_pipeline(muset_opt);
        } else {
            // use an input text matrix or a previous kmtricks directory
            spdlog::info(fmt::format("Using input k-mer matrix: {}", (muset_opt->in_matrix).c_str()));
            bool is_txt_input = fs::is_regular_file(muset_opt->in_matrix);
            // text matrix
            if(is_txt_input) {
                std::string kmer; kmat::TextMatrixReader mat(muset_opt->in_matrix);
                if (!mat.read_kmer(kmer)) {
                    throw std::runtime_error("Empty input matrix");
                }
                muset_opt->kmer_size = kmer.size();
                spdlog::debug(fmt::format("input matrix k-mer size: {}", muset_opt->kmer_size));
                muset_opt->kmer_matrix = muset_opt->in_matrix;
            }
            // kmtricks directory
            else if(!is_txt_input && kmat::is_kmtricks_dir(muset_opt->in_matrix)) {
                km::KmDir::get().init(muset_opt->in_matrix, "", false);
                Storage* config_storage = StorageFactory(STORAGE_FILE).load(km::KmDir::get().m_config_storage);
                LOCAL(config_storage);
                Configuration config = Configuration();
                config.load(config_storage->getGroup("gatb"));
                muset_opt->kmer_size = config._kmerSize;
                muset_opt->kmer_matrix = muset_opt->in_matrix;
            }
            else {
                throw std::runtime_error(fmt::format("input is neither a text file nor a valid kmtricks directory"));
            }
        }

        spdlog::info(fmt::format("Filtering k-mer matrix"));
        muset_opt->filtered_matrix = muset_opt->out_dir/"matrix.filtered.mat";
        kmat_filter(muset_opt);

        spdlog::info(fmt::format("Writing k-mers in FASTA format"));
        muset_opt->filtered_kmers = muset_opt->out_dir/"matrix.filtered.fasta";
        kmat_fasta(muset_opt);

        if(fs::is_empty(muset_opt->filtered_kmers)) {
            muset_opt->remove_temp_files();
            throw std::runtime_error("Filtered k-mer matrix is empty (filters were probably too strict).");
        }

        spdlog::info(fmt::format("Building unitigs"));
        muset_opt->unitigs = muset_opt->out_dir/"unitigs";
        ggcat(muset_opt);

        spdlog::info(fmt::format("Filtering unitigs"));
        muset_opt->filtered_unitigs = muset_opt->out_dir/"unitigs.fa";
        kmat_fafmt(muset_opt);

        if(fs::is_empty(muset_opt->filtered_unitigs)) {
            muset_opt->remove_temp_files();
            throw std::runtime_error("No unitig retained to build the output matrix (filters were probably too strict).");
        }

        spdlog::info(fmt::format("Building unitig matrix"));
        muset_opt->unitig_prefix = muset_opt->out_dir/"unitigs";
        kmat_unitig(muset_opt);

        spdlog::debug(fmt::format("Removing temporary files"));
        muset_opt->remove_temp_files();
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
