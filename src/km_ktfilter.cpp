#include <algorithm>
#include <filesystem>
#include <memory>
#include <numeric>
#include <string>

#include <fmt/format.h>

#define WITH_KM_IO
#include <kmtricks/public.hpp>

#include "common.h"

namespace fs = std::filesystem;

struct filter_options {
  fs::path matrices_dir;
  fs::path filtered_dir;
  fs::path output;
  std::size_t min_zeros{10};
  std::size_t min_nz{10};
  std::size_t min_abund{1};
  double min_zero_frac{0.1};
  double min_nz_frac{0.1};
  bool min_zero_frac_set{false};
  bool min_nz_frac_set{false};
  uint32_t kmer_size{31};
  std::size_t nb_threads{1};
};

void print_filter_usage() {
  filter_options opt;
  fmt::print("Usage: kmat_tools ktfilter [options] <kmtricks_run_dir>\n\n");
  fmt::print("Filter a kmtricks matrix by selecting k-mers that are potentially differential.\n\n");
  fmt::print("Options:\n");
  fmt::print("  -a INT    min abundance to define a k-mer as present in a sample [{}]\n", opt.min_abund);
  fmt::print("  -n INT    min number of samples for which a k-mer should be absent [{}]\n", opt.min_zeros);
  fmt::print("  -f FLOAT  fraction of samples for which a k-mer should be absent (overrides -n)\n");
  fmt::print("  -N INT    min number of samples for which a k-mer should be present [{}]\n", opt.min_nz);
  fmt::print("  -F FLOAT  fraction of samples for which a k-mer should be present (overrides -N)\n");
  fmt::print("  -o FILE   output filtered matrix to FILE [stdout]\n");
  fmt::print("  -t INT    number of threads [stdout]\n");
  fmt::print("  -w PATH   working directory for temporary files [\"<kmtricks_run_dir>/matrices_filtered\"]\n");
  fmt::print("  -h        print this help message\n");
}

template<size_t MAX_K>
class FilterTask : public km::ITask
{
  using count_type = typename km::selectC<DMAX_C>::type;

public:
  FilterTask(std::string &input, std::string &output, std::size_t &nb_kmers, std::size_t &nb_retained, filter_options &opts, bool compress = true)
    : km::ITask(4, false), m_input(input), m_output(output), m_nb_kmers(nb_kmers), m_nb_retained(nb_retained), m_opts(opts), m_compress(compress)
  {}

  void preprocess() {}
  void postprocess() {}

  void exec()
  {
    km::MatrixReader reader(m_input);
    km::Kmer<MAX_K> kmer; kmer.set_k(m_opts.kmer_size);
    
    std::size_t n_samples{reader.infos().nb_counts};
    std::vector<count_type> counts(n_samples);
    
    km::MatrixWriter<8192> writer(m_output,
      m_opts.kmer_size,
      reader.infos().count_slots,
      n_samples,
      reader.infos().id,
      reader.infos().partition,
      m_compress);

    while (reader.template read<MAX_K, DMAX_C>(kmer, counts)) {
      m_nb_kmers++;
      std::size_t n_zeros{0};
      std::size_t n_present{0};
      for (auto& c : counts) {
        if(c >= m_opts.min_abund){ 
          n_present++;
        } else { 
          n_zeros++; 
        }
      }

      bool enough_zeros = (m_opts.min_zero_frac_set && n_zeros >= m_opts.min_zero_frac*n_samples) || (!m_opts.min_zero_frac_set && n_zeros >= m_opts.min_zeros);
      bool enough_nz = (m_opts.min_nz_frac_set && n_present >= m_opts.min_nz_frac*n_samples) || (!m_opts.min_nz_frac_set && n_present >= m_opts.min_nz);
      if(enough_zeros && enough_nz) {
        m_nb_retained++;
        writer.template write<MAX_K, DMAX_C>(kmer,counts);
      }
    }
  }

private:
  std::string& m_input;
  std::string& m_output;
  std::size_t& m_nb_kmers;
  std::size_t& m_nb_retained;
  filter_options& m_opts;
  bool m_compress;
};

template<size_t MAX_K>
struct filter_functor {

  using count_type = typename km::selectC<DMAX_C>::type;

  void operator()(filter_options &opts)
  {
    std::vector<std::string> matrix_paths;
    std::vector<std::string> filtered_paths;
    for (auto const& entry : std::filesystem::directory_iterator{opts.matrices_dir}) {
      if(fs::is_regular_file(entry)) {
        matrix_paths.push_back(entry.path());
        filtered_paths.push_back(opts.filtered_dir/entry.path().filename());
      }
    }

    std::size_t nb_threads = std::min(opts.nb_threads,matrix_paths.size());
    km::TaskPool pool(nb_threads);
    std::vector<std::size_t> nb_total_kmers(matrix_paths.size(),0);
    std::vector<std::size_t> nb_retained(matrix_paths.size(),0);
    for (std::size_t i=0; i < matrix_paths.size(); i++) {
      pool.add_task(std::make_shared<FilterTask<MAX_K>>(matrix_paths[i], filtered_paths[i], nb_total_kmers[i], nb_retained[i], opts));
    }
    pool.join_all();

    km::MatrixFileAggregator<MAX_K,DMAX_C> mfa(filtered_paths, opts.kmer_size);
    opts.output.empty() ? mfa.write_as_text(std::cout) : mfa.write_as_text(opts.output);

    fmt::print(stderr, "[info] {} total kmers\n", std::reduce(nb_total_kmers.begin(), nb_total_kmers.end()));
    fmt::print(stderr, "[info] {} retained kmers\n", std::reduce(nb_retained.begin(), nb_retained.end()));
  }
};


int main_ktfilter(int argc, char **argv) {

  filter_options opts;

  int c;
  while ((c = getopt(argc, argv, "a:f:F:n:N:o:t:w:h")) != -1) {
    switch (c) {
      case 'a':
        opts.min_abund = strtoul(optarg, NULL, 10);
        break;
      case 'n':
        opts.min_zeros = strtoul(optarg, NULL, 10);
        break;
      case 'N':
        opts.min_nz = strtoul(optarg, NULL, 10);
        break;
      case 'o':
        opts.output = optarg;
        break;
      case 'f':
        opts.min_zero_frac_set = true;
        opts.min_zero_frac = atof(optarg);
        break;
      case 'F':
        opts.min_nz_frac_set = true;
        opts.min_nz_frac = atof(optarg);
        break;
      case 'h':
        print_filter_usage();
        return 0;
      case 't':
        opts.nb_threads = std::max(1, std::stoi(optarg));
        break;
      case 'w':
        opts.filtered_dir = optarg;
        break;
      case '?':
        return 1;
      default:
        abort();
    }
  }

  if(argc-optind != 1) {
    print_filter_usage();
    return 0;
  }

  const fs::path kmtricks_dir{argv[optind]};
  opts.matrices_dir = kmtricks_dir/"matrices";

  if(!fs::is_directory(kmtricks_dir)) {
    fmt::print(stderr, "[error] input directory \"{}\" does not exist\n", kmtricks_dir.c_str());
    return 1;
  }

  if(!fs::is_directory(opts.matrices_dir)) {
    fmt::print(stderr, "[error] kmtricks \"matrices\" directory does not exist\n");
    return 1;
  }

  // empty kmtricks matrix
  if (fs::is_empty(opts.matrices_dir)) {
    fmt::print(stderr, "[error] kmtricks matrices directory is empty.");
    return 1;
  }

  // retrieve k-mer size from one of the matrix files
  for (auto const& entry : std::filesystem::directory_iterator{opts.matrices_dir}) {
    km::MatrixReader reader(entry.path());
    opts.kmer_size = reader.infos().kmer_size;
    break;
  }
  
  if (opts.filtered_dir.empty()) { opts.filtered_dir = kmtricks_dir/"matrices_filtered"; }

  // create working directory and check it is not the same as the matrices directory
  fs::create_directories(opts.filtered_dir);
  if (fs::equivalent(opts.matrices_dir,opts.filtered_dir)) {
    fmt::print(stderr, "[error] working directory cannot be the same as the matrices directory");
    std::exit(EXIT_FAILURE);
  }
  
  try 
  { 
    km::const_loop_executor<0, KMER_N>::exec<filter_functor>(opts.kmer_size, opts);
  }
  catch (const km::km_exception &e)
  {
    fmt::print(stderr, "[exception] {} - {}", e.get_name(), e.get_msg());
    std::exit(EXIT_FAILURE);
  }

  return 0;
}
