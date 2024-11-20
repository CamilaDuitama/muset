#pragma once

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <kseq++/seqio.hpp>
#include "../external/sshash/dictionary.hpp"

#include <kmat_tools/cli/unitig.h>
#include <kmat_tools/matrix.h>
#include <kmat_tools/task.h>
#include <kmat_tools/utils.h>

namespace fs = std::filesystem;


namespace kmat {

int main_unitig(unitig_opt_t opt)
{
    // input validation

    fs::path unitig_path = opt->inputs[0];
    if(!fs::is_regular_file(unitig_path)) {
        throw std::runtime_error(fmt::format("unitig file \"{}\" does not exist", unitig_path.c_str()));
    }

    fs::path matrix_path = opt->inputs[1];
    if(!fs::is_regular_file(matrix_path)) {
        throw std::runtime_error(fmt::format("k-mer matrix file \"{}\" does not exist", matrix_path.c_str()));
    }

    if(opt->mini_size >= opt->kmer_size) {
        throw std::runtime_error("minimizer size must be smaller than k-mer size");
    }

    spdlog::info("building k-mer dictionary");

    sshash::dictionary kmer_dict;
    {
        // std::ofstream ofs("sshash.log", std::ios::out);
        // std::streambuf *coutbuf = std::cout.rdbuf();
        // if (ofs.good()) {
        //   std::cout.rdbuf(ofs.rdbuf());
        // }

        sshash::build_configuration build_config;
        build_config.k = opt->kmer_size;
        build_config.m = opt->mini_size;
        build_config.c = 5.0;
        build_config.pthash_threads = opt->nb_threads;
        build_config.canonical_parsing = true;
        build_config.verbose = false;

        kmer_dict.build(unitig_path, build_config);

        // std::cout.rdbuf(coutbuf);
    }

    spdlog::debug(fmt::format("k-mer processed: {}", kmer_dict.size()));
    spdlog::debug(fmt::format("unitigs processed: {}", kmer_dict.num_contigs()));

    spdlog::info("aggregating k-mer counts");

    TextMatrixReader mat(matrix_path);

    std::string kmer;
    std::vector<uint32_t> kmer_counts;
    bool has_kmer = mat.read_kmer_counts(kmer,kmer_counts);
    
    std::size_t nb_samples {has_kmer ? kmer_counts.size() : 0};
    spdlog::debug(fmt::format("samples: {}", nb_samples));

    std::vector<std::vector<std::pair<uint32_t,uint32_t>>> unitig_counts(kmer_dict.num_contigs());
    for(auto& samples: unitig_counts) { samples.resize(nb_samples); }

    while(has_kmer) {

        auto res = kmer_dict.lookup_advanced(kmer.c_str());
        if (res.kmer_id == sshash::constants::invalid_uint64) {
            has_kmer = mat.read_kmer_counts(kmer,kmer_counts);
            continue;
        }

        std::size_t utg_id = res.contig_id;
        auto& samples = unitig_counts[utg_id];
        
        size_t idx{0};
        for(auto& [nb_present,abundance_sum]: samples) {
            uint32_t num = kmer_counts[idx];
            nb_present = add_sat(nb_present, uint32_t{num > 0});
            abundance_sum = add_sat(abundance_sum, num);
            idx++;
        }

        has_kmer = mat.read_kmer_counts(kmer,kmer_counts);
    }

    spdlog::info("writing unitig matrix");

    std::string matrix_file = fmt::format("{}.abundance.mat", opt->prefix);
    std::string frac_file = fmt::format("{}.frac.mat", opt->prefix);

    std::ofstream mat_ofs(matrix_file);
    if (!mat_ofs.good()) { throw std::runtime_error(fmt::format("cannot open output file \"{}\" for writing", matrix_file)); }
    mat_ofs << std::fixed << std::setprecision(2);

    std::ofstream frac_ofs;
    if (opt->write_frac_matrix) {
        frac_ofs.open(frac_file);
        if (!frac_ofs.good()) { throw std::runtime_error(fmt::format("cannot open output file \"{}\" for writing", frac_file)); }
    }

    klibpp::KSeq unitig;
    klibpp::SeqStreamIn utg_ssi(unitig_path.c_str());
    for(uint64_t utg_id=0; utg_ssi >> unitig; utg_id++) {  

        mat_ofs << (opt->write_seq ? unitig.seq : unitig.name);
        auto& samples = unitig_counts[utg_id];
        
        std::size_t utg_nb_kmers {unitig.seq.length() - opt->kmer_size + 1};
        for(auto [nb_present,abundance_sum] : samples) {
            double frac = (1.0 * nb_present)/utg_nb_kmers;
            double avg_abundance = (1.0 * abundance_sum)/utg_nb_kmers;
            mat_ofs << " " << (frac >= opt->min_frac ? avg_abundance : 0.0);
            if(opt->write_frac_matrix) { frac_ofs << " " << frac; }
        }

        mat_ofs << "\n";
        if(opt->write_frac_matrix) { frac_ofs << "\n"; }
    }

    mat_ofs.close();
    spdlog::info(fmt::format("abundance unitig matrix written to \"{}\"", matrix_file));

    if(opt->write_frac_matrix) {
        frac_ofs.close();
        spdlog::info(fmt::format("fraction unitig matrix written to \"{}\"", frac_file));
    }

    return 0;
}

};