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

#include <kmat_tools/cli.h>
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
        std::string sshash_logfile = fmt::format("{}.sshash.log", opt->prefix);
        std::ofstream ofs(sshash_logfile, std::ios::out);
        std::streambuf *coutbuf = std::cout.rdbuf();
        if (ofs.good()) {
          std::cout.rdbuf(ofs.rdbuf());
        }

        sshash::build_configuration build_config;
        build_config.k = opt->kmer_size;
        build_config.m = opt->mini_size;
        build_config.c = 5.0;
        build_config.pthash_threads = opt->nb_threads;
        build_config.canonical_parsing = true;
        build_config.verbose = false;
        kmer_dict.build(unitig_path, build_config);

        std::cout.rdbuf(coutbuf);
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
        frac_ofs << std::fixed << std::setprecision(2);
    }

    klibpp::KSeq unitig;
    klibpp::SeqStreamIn utg_ssi(unitig_path.c_str());
    for(uint64_t utg_id=0; utg_ssi >> unitig; utg_id++) {  

        mat_ofs << (opt->write_seq ? unitig.seq : unitig.name);
        if(opt->write_frac_matrix) { frac_ofs << (opt->write_seq ? unitig.seq : unitig.name); }
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


kmat_opt_t unitig_cli(std::shared_ptr<bc::Parser<1>> cli, unitig_opt_t opt)
{
    bc::cmd_t unitig = cli->add_command("unitig", "Create a unitig matrix.");

    unitig->add_group("main options", "");

    unitig->add_param("-k/--kmer-size", fmt::format("k-mer size [8,{}].", KL[MUSET_KMER_N-1]-1))
        ->meta("INT")
        ->def("31")
        ->checker(bc::check::f::range(8, KL[MUSET_KMER_N-1]-1))
        ->setter(opt->kmer_size);

    unitig->add_param("-p/--prefix", "output files prefix.")
        ->meta("FILE")
        ->def("out")
        ->setter(opt->prefix);

    unitig->add_param("-f/--min-frac", "set average abundance to 0 if k-mer fraction is below this threshold [0,1].")
        ->meta("FLOAT")
        ->def("0.0")
        ->checker(bc::check::f::range(0.0, 1.0))
        ->setter(opt->min_frac);

    unitig->add_param("--out-frac", "output an additional matrix containing k-mer fractions.")
        ->as_flag()
        ->setter(opt->write_frac_matrix);

    unitig->add_param("-s/--write-seq", "write the unitig sequence instead of the identifier in the output matrix")
        ->as_flag()
        ->setter(opt->write_seq);

    unitig->add_group("other options", "");

    unitig->add_param("-m/--minimizer-size", "minimizer size")
        ->meta("INT")
        ->def("15")
        ->setter(opt->mini_size);

    unitig->add_param("-t/--threads", "number of threads.")
        ->meta("INT")
        ->def("4")
        ->checker(bc::check::is_number)
        ->setter(opt->nb_threads);

    unitig->add_param("-h/--help", "show this message and exit.")
         ->as_flag()
         ->action(bc::Action::ShowHelp);
    
    unitig->add_param("-v/--version", "show version and exit.")
         ->as_flag()
         ->action(bc::Action::ShowVersion);

    unitig->set_positionals(2, "<unitigs.fasta> <kmer_matrix>", "a unitig fasta file and a text-based k-mer matrix");

    return opt;
}

};