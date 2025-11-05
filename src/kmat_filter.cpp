#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#define WITH_KM_IO
#include <kmtricks/public.hpp>

#include <kmat_tools/cmd/filter.h>
#include <kmat_tools/matrix.h>
#include <kmat_tools/task.h>
#include <kmat_tools/utils.h>

namespace fs = std::filesystem;


namespace kmat {

bool should_skip_filter(filter_opt_t opt) {
    // Check if filtering would let everything through
    bool no_absent_filter = (!opt->min_nb_absent_set && opt->min_frac_absent == 0.0)
                         || (opt->min_nb_absent_set && opt->min_nb_absent == 0);
    
    bool no_present_filter = (!opt->min_nb_present_set && opt->min_frac_present == 0.0)
                          || (opt->min_nb_present_set && opt->min_nb_present == 0);
    
    return (no_absent_filter && no_present_filter);
}

int kmat_basic_filter(fs::path input, filter_opt_t opt) {

    // Optimization: if no filtering is needed, just copy the file
    if (should_skip_filter(opt) && !(opt->output).empty()) {
        spdlog::info(fmt::format("No filtering needed - copying matrix"));
        std::ifstream src(input, std::ios::binary);
        std::ofstream dst(opt->output, std::ios::binary);
        dst << src.rdbuf();
        
        // Count lines for reporting
        src.clear();
        src.seekg(0);
        size_t nb_kmers = 0;
        std::string line;
        while (std::getline(src, line)) { nb_kmers++; }
        
        spdlog::info(fmt::format("{}/{} k-mers retained (no filtering)", nb_kmers, nb_kmers));
        return 0;
    }

    TextMatrixReader reader(input);

    std::ostream* fpout = &std::cout;
    std::ofstream ofs;
    if(!(opt->output).empty()) {
        ofs.open((opt->output).c_str());
        if(!ofs.good()) { throw std::runtime_error(fmt::format("cannot open output file {}", (opt->output).c_str())); }
        fpout = &ofs;
    }

    size_t nb_samples{0};
    size_t nb_kmers{0};
    size_t nb_retained{0};

    std::string kmer;
    std::vector<size_t> counts;

    while(reader.read_kmer_counts(kmer,counts)) {

        nb_kmers++;
        if (nb_kmers == 1) { nb_samples = counts.size(); }

        if (nb_samples != counts.size()) {
            throw std::runtime_error(fmt::format("inconsistent number of samples at line {}: found {}, expected {}", reader.line_count(), counts.size(), nb_samples));
        }

        size_t nb_absent{0};
        size_t nb_present{0};
        for (auto value: counts) {
            if(value >= opt->min_abundance){
                nb_present++;
            } else {
                nb_absent++;
            }
        }

        bool enough_absent = (!opt->min_nb_absent_set && nb_absent >= opt->min_frac_absent * nb_samples)
            || (opt->min_nb_absent_set && nb_absent >= opt->min_nb_absent);    
        bool enough_present = (!opt->min_nb_present_set && nb_present >= opt->min_frac_present * nb_samples)
            || (opt->min_nb_present_set && nb_present >= opt->min_nb_present);
        if (enough_absent && enough_present) {
            nb_retained++;
            *fpout << reader.line() << "\n";
        }
    }

    spdlog::info(fmt::format("{} samples", nb_samples));
    spdlog::info(fmt::format("{}/{} k-mers retained", nb_retained, nb_kmers));

    if(!(opt->output).empty()) {
        ofs.close();
    }

    return 0;
}


template<size_t MAX_K>
struct kmtricks_matrix_filter {
    
    using count_type = typename km::selectC<DMAX_C>::type;

    void operator()(filter_opt_t opts)
    {
        std::vector<std::string> matrix_paths;
        std::vector<std::string> filtered_paths;
        for (auto const& entry : fs::directory_iterator{opts->matrices_dir}) {
            if(fs::is_regular_file(entry)) {
                matrix_paths.push_back(entry.path());
                filtered_paths.push_back((opts->filtered_dir)/entry.path().filename());
            }
        }

        // Check if we can skip filtering
        if (should_skip_filter(opts)) {
            spdlog::info(fmt::format("No filtering needed - aggregating matrices"));
            // Directly aggregate without filtering
            km::MatrixFileAggregator<MAX_K,DMAX_C> mfa(matrix_paths, opts->kmer_size);
            (opts->output).empty() ? mfa.write_as_text(std::cout) : mfa.write_as_text(opts->output);
            spdlog::info(fmt::format("All k-mers retained (no filtering)"));
        } else {
            std::size_t nb_threads = std::min(opts->nb_threads, matrix_paths.size());
            km::TaskPool pool(nb_threads);
            std::vector<std::size_t> nb_total_kmers(matrix_paths.size(),0);
            std::vector<std::size_t> nb_retained(matrix_paths.size(),0);
            for (std::size_t i=0; i < matrix_paths.size(); i++) {
                pool.add_task(std::make_shared<FilterTask<MAX_K>>(matrix_paths[i], filtered_paths[i], nb_total_kmers[i], nb_retained[i], opts));
            }
            pool.join_all();

            km::MatrixFileAggregator<MAX_K,DMAX_C> mfa(filtered_paths, opts->kmer_size);
            (opts->output).empty() ? mfa.write_as_text(std::cout) : mfa.write_as_text(opts->output);

            auto retained_kmers = std::reduce(nb_retained.begin(), nb_retained.end());
            auto total_kmers = std::reduce(nb_total_kmers.begin(), nb_total_kmers.end());
            spdlog::info(fmt::format("{}/{} kmers retained", retained_kmers, total_kmers));
        }
    }
};


int main_filter(filter_opt_t opt)
{
    // input validation

    fs::path input = opt->inputs[0];

    bool is_txt_input = fs::is_regular_file(input);
    bool is_dir_input = !is_txt_input && is_kmtricks_dir(input);

    if(!is_txt_input && !is_dir_input) {
        throw std::runtime_error(fmt::format("input is neither a file nor a valid kmtricks directory"));
    }

    if (is_txt_input) {
        spdlog::info(fmt::format("filtering text matrix: {}", input.c_str()));
        return kmat_basic_filter(input, opt);
    }

    // here the input should be a kmtricks directory

    const fs::path kmtricks_dir{input};
    opt->matrices_dir = kmtricks_dir/"matrices";
    opt->filtered_dir = kmtricks_dir/"matrices_filtered";

    spdlog::info(fmt::format("filtering kmtricks matrix: {}", kmtricks_dir.c_str()));

    // empty kmtricks matrix
    if (fs::is_empty(opt->matrices_dir)) {
        throw std::runtime_error(fmt::format("kmtricks matrices directory is empty"));
    }

    // retrieve k-mer size from one of the matrix files
    for (auto const& entry : fs::directory_iterator{opt->matrices_dir}) {
        km::MatrixReader reader(entry.path());
        opt->kmer_size = reader.infos().kmer_size;
        break;
    }

    fs::create_directories(opt->filtered_dir);
    km::const_loop_executor<0, KMER_N>::exec<kmtricks_matrix_filter>(opt->kmer_size, opt);
    if (!opt->keep_tmp) { fs::remove_all(opt->filtered_dir); }

    return 0;
}

};