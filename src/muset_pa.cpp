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
#include <kmtricks/loop_executor.hpp>

#include <kmat_tools/cmd/convert.h>
#include <kmat_tools/cmd/fafmt.h>

#include "muset_pa_cli.h"


void init_logger(const fs::path &dir) {
    auto cerr_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
    cerr_sink->set_pattern("[%Y-%m-%d %H:%M:%S] [%^%l%$] %v");

    auto now = std::chrono::system_clock::now();
    auto local_time = std::chrono::system_clock::to_time_t(now);
    std::stringstream ss; ss << std::put_time(std::localtime(&local_time), "%Y%m%d_%H%M%S");
    auto muset_log = fmt::format("{}/muset_pa_{}.log", dir.c_str(), ss.str());
    auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(muset_log, true);
    file_sink->set_pattern("[%Y-%m-%d %H:%M:%S] [%^%l%$] %v");

    auto combined_logger = std::make_shared<spdlog::logger>("muset", spdlog::sinks_init_list({cerr_sink,file_sink}));
    spdlog::set_default_logger(combined_logger);
}

void print_options(muset::muset_pa_options_t opt) {

    spdlog::info(fmt::format("input file (--file): {}", (opt->fof).c_str()));
    spdlog::info(fmt::format("output directory (-o): {}", (opt->out_dir).c_str()));
    
    spdlog::info(fmt::format("k-mer size (-k): {}", opt->kmer_size));
    spdlog::info(fmt::format("minimum abundance (-a): {}", opt->min_abundance));
    spdlog::info(fmt::format("minimum unitig length (-l): {}", opt->min_utg_len));
    spdlog::info(fmt::format("minimum unitig fraction set: {}", opt->min_utg_frac_set));
    if(opt->min_utg_frac_set) { spdlog::info(fmt::format("minimum unitig fraction (-r): {}", opt->min_utg_frac)); }
    spdlog::info(fmt::format("write unitig sequence (-s): {}", opt->write_utg_seq));
    spdlog::info(fmt::format("minimizer size (-m): {}", opt->mini_size));

    spdlog::info(fmt::format("keep temporary files (--keep-temp): {}", opt->keep_tmp));
    spdlog::info(fmt::format("threads (-t): {}", opt->nb_threads));
}

void ggcat_build(muset::muset_pa_options_t opt) {
    std::string ggcat_tempdir = opt->out_dir/"ggcat_build";
    std::string ggcat_keeptmp = opt->keep_tmp ? "--keep-temp-files" : "";
    std::string ggcat_logfile = opt->out_dir/"ggcat_build.log";
    std::string ggcat_cmd = fmt::format("ggcat build -k {} -s {} --minimizer-length {} --colors -j {} -l {} -o {} {} --temp-dir {} &> {}",
        opt->kmer_size, // -k
        opt->min_abundance, // -s
        opt->mini_size, // --minimizer-length
        opt->nb_threads, // -j
        (opt->fof).c_str(), // -l
        (opt->unitigs).c_str(), // -o
        ggcat_keeptmp, // --keep-temp-files ?
        ggcat_tempdir, // --temp-dir
        ggcat_logfile);

    spdlog::debug(fmt::format("Running command: {}", ggcat_cmd));    
    auto ret = std::system(ggcat_cmd.c_str());
    if(ret != 0) {
        throw std::runtime_error(fmt::format("Command failed: {}\nSee log at {}", ggcat_cmd, ggcat_logfile));
    }
}

void ggcat_dump_colors(muset::muset_pa_options_t opt) {
    std::string colormap = fmt::format("{}.colors.dat", (opt->unitigs).c_str());
    std::string ggcat_logfile = opt->out_dir/"ggcat_dump-colors.log";
    std::string ggcat_cmd = fmt::format("ggcat dump-colors {} {} &> {}",
        colormap, (opt->colors_json).c_str(), ggcat_logfile);

    spdlog::debug(fmt::format("Running command: {}", ggcat_cmd));    
    auto ret = std::system(ggcat_cmd.c_str());
    if(ret != 0) {
        throw std::runtime_error(fmt::format("Command failed: {}\nSee log at {}", ggcat_cmd, ggcat_logfile));
    }
}

void ggcat_query(muset::muset_pa_options_t opt) {
    std::string ggcat_tempdir = opt->out_dir/"ggcat_query";
    std::string ggcat_keeptmp = opt->keep_tmp ? "--keep-temp-files" : "";
    std::string ggcat_logfile = opt->out_dir/"ggcat_query.log";
    std::string ggcat_cmd = fmt::format("ggcat query --colors -k {} --minimizer-length {} -j {} -o {} {} --temp-dir {} {} {} &> {}",
        opt->kmer_size, // -k
        opt->mini_size, // --minimizer-length
        opt->nb_threads, // -j
        (opt->query_json).c_str(), // -o
        ggcat_keeptmp, // --keep-temp-files ?
        ggcat_tempdir, // --temp-dir
        (opt->filtered_unitigs).c_str(),
        (opt->filtered_unitigs).c_str(),
        ggcat_logfile);

    spdlog::debug(fmt::format("Running command: {}", ggcat_cmd));    
    auto ret = std::system(ggcat_cmd.c_str());
    if(ret != 0) {
        throw std::runtime_error(fmt::format("Command failed: {}\nSee log at {}", ggcat_cmd, ggcat_logfile));
    }
}

void kmat_fafmt(muset::muset_pa_options_t opt) {

    auto fafmt_opt = std::make_shared<kmat::fafmt_options>();
    fafmt_opt->output = opt->filtered_unitigs;
    fafmt_opt->min_length = opt->min_utg_len;
    (fafmt_opt->inputs).push_back(opt->unitigs);

    kmat::main_fafmt(fafmt_opt);
}

void kmat_convert(muset::muset_pa_options_t opt) {

    auto convert_opt = std::make_shared<kmat::convert_options>();

    convert_opt->ap_flag = opt->min_utg_frac_set;
    convert_opt->min_frac_set = opt->min_utg_frac_set;
    convert_opt->min_frac = opt->min_utg_frac;

    convert_opt->out_write_seq = opt->write_utg_seq;
    // convert_opt->no_header = true;
    convert_opt->no_header = true;
    convert_opt->out_csv = false;
    convert_opt->out_fname = opt->unitig_matrix;

    (convert_opt->inputs).push_back(opt->filtered_unitigs);
    (convert_opt->inputs).push_back(opt->colors_json);
    (convert_opt->inputs).push_back(opt->query_json);

    kmat::main_convert(convert_opt);
}

int main(int argc, char* argv[])
{
    muset::musetPaCli cli("muset_pa", "a pipeline for building a presence-absence unitig matrix from a list of FASTA/FASTQ files.", PROJECT_VER, "");
    auto opt = cli.parse(argc, argv);

    try
    {
        // check ggcat dependency
        if (std::system("ggcat --version &>/dev/null") != 0) {
            throw std::runtime_error("ggcat: command not found");
        }

        // check parameters consistency
        opt->sanity_check();
        
        // minimum unitig length is set to 2k-1 if no threshold provided
        if (!opt->min_utg_len_set) {
            opt->min_utg_len = 2 * opt->kmer_size - 1;
        }

        // create main output directory
        fs::create_directories(opt->out_dir);

        // initialize the logger just before running the pipeline
        init_logger(opt->out_dir);

        // print values of main options
        spdlog::info(fmt::format("Running muset_pa {}", PROJECT_VER));
        spdlog::info("---------------");
        print_options(opt);
        spdlog::info("---------------");

        // muset_pa pipeline

        spdlog::info(fmt::format("Building unitigs"));
        opt->unitigs = opt->out_dir/"unitigs";
        ggcat_build(opt);

        opt->colors_json = opt->out_dir/"unitigs.jsonl";
        ggcat_dump_colors(opt);

        spdlog::info(fmt::format("Filtering unitigs"));
        opt->filtered_unitigs = opt->out_dir/"unitigs.fa";
        kmat_fafmt(opt);

        if(fs::is_empty(opt->filtered_unitigs)) {
            throw std::runtime_error("No unitig retained to build the output matrix (filters were probably too strict).");
        }

        opt->query_json = opt->out_dir/"unitigs.query.jsonl";
        ggcat_query(opt);

        spdlog::info(fmt::format("Building unitig matrix"));
        opt->unitig_matrix = opt->out_dir/"unitigs.mat";
        kmat_convert(opt);
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
