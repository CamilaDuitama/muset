#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <sstream>

#include <fmt/format.h>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <kmtricks/cmd.hpp>
#include <kmtricks/loop_executor.hpp>

#include <kmat_tools/cmd/convert.h>
#include <kmat_tools/cmd/fafmt.h>

#include "muset_pa_cli.h"


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

    // spdlog::set_level(spdlog::level::debug);
    auto cerr_logger = spdlog::stderr_color_mt("muset");
    cerr_logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");
    spdlog::set_default_logger(cerr_logger);

    try
    {
        // check ggcat dependency
        if (std::system("ggcat --version &>/dev/null") != 0) {
            throw std::runtime_error("ggcat: command not found");
        }

        auto opt = cli.parse(argc, argv);
        opt->sanity_check();
        
        // minimum unitig length is set to 2k-1 if no threshold provided
        if (!opt->min_utg_len_set) {
            opt->min_utg_len = 2 * opt->kmer_size - 1;
        }

        fs::create_directories(opt->out_dir);

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
