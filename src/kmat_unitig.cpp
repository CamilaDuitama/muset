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

#include <kmat_tools/cmd/unitig.h>
#include <kmat_tools/matrix.h>
#include <kmat_tools/task.h>
#include <kmat_tools/utils.h>

#include <kmat_tools/aggregator.h>
#include <kmat_tools/matrix_writer.h>

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

    size_t number_unitigs {kmer_dict.num_contigs()};

    std::unique_ptr<Aggregator> aggregator;
    if (opt->abundance_metric == "mean") {
        spdlog::info(fmt::format("Computing mean ({}) for unitigs.", opt->abundance_metric));
        aggregator = std::make_unique<MeanAggregator>(nb_samples, number_unitigs, opt->min_frac );
    }
    else { // median
        spdlog::info(fmt::format("Computing median ({}) for unitigs.", opt->abundance_metric));
        aggregator = std::make_unique<MedianAggregator>(nb_samples, number_unitigs, opt->min_frac );
    }

    while(has_kmer) {
        auto res = kmer_dict.lookup_advanced(kmer.c_str());
        if (res.kmer_id == sshash::constants::invalid_uint64) {
            has_kmer = mat.read_kmer_counts(kmer,kmer_counts);
            continue;
        }

        std::size_t utg_id = res.contig_id;
        aggregator->process_kmer(utg_id, kmer_counts);

        has_kmer = mat.read_kmer_counts(kmer,kmer_counts);
    }

    spdlog::info("writing unitig matrix");

    std::unique_ptr<MatrixWriter> writer;
    if (opt->output_format == "txt") {
        // NORMAL TEXT FILE
        writer = std::make_unique<TextMatrixWriter>(opt->prefix, opt->write_frac_matrix );
    } else if (opt->output_format == "tsv") {
        // TSV COMPRESSED FOR DOWNSTREAM IN MEMORY
        writer = std::make_unique<CompressedTSVMatrixWriter>(opt->prefix, opt->write_frac_matrix, nb_samples);
    } else {
        // IF NOT RECOGNIZED DEFAULT TO TXT FOR BACKWARD COMPATIBILITY AND AVOID DISRUPTION
        spdlog::info(fmt::format("OUTPUT FORMAT {} NOT RECOGNIZED. DEFAULT TO TXT.", opt->output_format));
        writer = std::make_unique<TextMatrixWriter>(opt->prefix, opt->write_frac_matrix );
    }

    klibpp::KSeq unitig;
    klibpp::SeqStreamIn utg_ssi(unitig_path.c_str());

    std::vector<double> utg_abundances(nb_samples);
    std::vector<double> utg_fractions(nb_samples);

    for(uint64_t utg_id=0; utg_ssi >> unitig; utg_id++) {
        std::size_t utg_nb_kmers {unitig.seq.length() - opt->kmer_size + 1};
        // FILL SAMPLES VECTOR WITH ABUNDANCE FRACTION FOR THE UTG
        for (size_t idx {0}; idx < nb_samples; idx++){
            auto [abundance, frac] = aggregator->get_abundance_fraction(utg_id, idx, utg_nb_kmers);
            utg_abundances[idx] = abundance;
            if (opt->write_frac_matrix) {utg_fractions[idx] = frac;}
        }
        // DUMP IT TO DISK
        std::string utg_identifier {opt->write_seq ? unitig.seq : unitig.name};
        writer->write_row(utg_identifier, utg_abundances, utg_fractions);
        // REPEAT
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

    unitig->add_param("-f/--min-frac", "set unitig average abundance to 0 if its k-mer fraction is below this threshold [0,1].")
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

    unitig->add_param("--abundance-metric", "metric to use for abundance: mean or median.")
    ->meta("STRING")
    ->def("mean") // Default to mean for backward compatibility
    ->checker(bc::check::f::in("mean|median"))
    ->setter(opt->abundance_metric);

    unitig->add_param("--output-format", "Output format can be either 'txt' or 'tsv' (tsv is gzip compressed).")
    ->meta("STRING")
    ->def("txt") // Default to txt for backward compatibility
    ->checker(bc::check::f::in("txt|tsv"))
    ->setter(opt->output_format);

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