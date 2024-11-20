#pragma once

#include <filesystem>
#include <memory>
#include <string>

#include <algorithm>
#include <bcli/bcli.hpp>

#include <kmtricks/exceptions.hpp>
#include <kmtricks/utils.hpp>

#define RECORD(ss, var) ss << #var << "=" << var << ", "

namespace fs = std::filesystem;


namespace muset {

using muset_options_t = std::shared_ptr<struct muset_options>;
using cli_t = std::shared_ptr<bc::Parser<0>>;



struct muset_options
{
  fs::path fof;
  fs::path in_matrix;
  fs::path out_dir;

  uint32_t kmer_size{0};
  uint32_t mini_size{0};
  uint32_t min_abundance{0};
  bool min_utg_len_set{false};
  uint32_t min_utg_len{0};
  bool write_utg_seq{false};
  
  bool no_kmer_filter{false};
  double min_frac_absent{0.1};
  double min_frac_present{0.1};
  bool min_nb_absent_set{false};
  int min_nb_absent{0};
  bool min_nb_present_set{false};
  int min_nb_present{0};
  
  bool keep_tmp{false};
  bool lz4{true};
  bool logan{false};

  int nb_threads{1};

  // intermediate files, not actual input parameters

  fs::path filtered_matrix;

  std::string display()
  {
    std::stringstream ss;
    
    RECORD(ss, fof);
    RECORD(ss, in_matrix);
    RECORD(ss, out_dir);

    RECORD(ss, kmer_size);
    RECORD(ss, mini_size);
    RECORD(ss, min_abundance);
    RECORD(ss, min_utg_len_set);
    RECORD(ss, min_utg_len);
    RECORD(ss, write_utg_seq);
    
    RECORD(ss, no_kmer_filter);
    RECORD(ss, min_frac_absent);
    RECORD(ss, min_frac_present);
    RECORD(ss, min_nb_absent_set);
    RECORD(ss, min_nb_absent);
    RECORD(ss, min_nb_present_set);
    RECORD(ss, min_nb_present);
    
    RECORD(ss, keep_tmp);
    RECORD(ss, lz4);
    RECORD(ss, logan);

    RECORD(ss, nb_threads);

    return ss.str();
  }

  void sanity_check()
  {

    // check input-file parameters
    if (in_matrix.empty() && fof.empty()) {
      throw std::runtime_error("either --file or --in-matrix should be provided");
    } else if (!in_matrix.empty() && !fof.empty()) {
      throw std::runtime_error("either --file or --in-matrix should be provided, but not both");
    }

    // check existance of the provided input file
    if (!fof.empty() && (!fs::is_regular_file(fof) || fs::is_empty(fof))) {
      throw std::runtime_error(fmt::format("input file \"{}\" does not exist", fof.c_str()));
    }
    if (!in_matrix.empty() && (!fs::is_regular_file(in_matrix) || fs::is_empty(in_matrix))) {
      throw std::runtime_error(fmt::format("input matrix \"{}\" is not a file or is empty", in_matrix.c_str()));
    }
    
    // --logan flag only valid for kmer-size equal to 31
    if ((logan) && (kmer_size != 31))
    {
      throw std::runtime_error("--logan available only for --kmer-size equal to 31");
    }

    if (mini_size >= kmer_size) {
      throw std::runtime_error("minimizer size must be smaller than k-mer size");
    }
  }

  void dump(const std::string& path)
  {
    std::ofstream out_opt(path, std::ios::out);
    km::check_fstream_good(path, out_opt);
    out_opt << display();
  }
};

class musetCli
{

public:
    musetCli(
        const std::string& name,
        const std::string& desc,
        const std::string& version,
        const std::string& authors
    );

    muset_options_t parse(int argc, char* argv[]);

private:
  cli_t cli{nullptr};
  muset_options_t muset_opt{nullptr};
};

muset_options_t muset_cli(std::shared_ptr<bc::Parser<0>> cli, muset_options_t options);

};  // namespace km
