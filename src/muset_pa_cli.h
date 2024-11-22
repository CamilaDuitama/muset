#pragma once

#include <filesystem>
#include <memory>
#include <string>

#include <algorithm>
#include <bcli/bcli.hpp>

namespace fs = std::filesystem;


namespace muset {

using muset_pa_cli_t = std::shared_ptr<bc::Parser<0>>;

struct muset_pa_options
{
    fs::path fof;
    fs::path out_dir;

    uint32_t kmer_size{0};
    uint32_t mini_size{0};
    uint32_t min_abundance{0};
    
    bool min_utg_len_set{false};
    uint32_t min_utg_len{0};

    bool min_utg_frac_set{false};
    double min_utg_frac{0.8};

    bool write_utg_seq{false};
    bool keep_tmp{false};
    int nb_threads{1};

    // intermediate files, not actual input parameters

    fs::path unitigs;
    fs::path filtered_unitigs;
    fs::path query_json;
    fs::path colors_json;
    fs::path unitig_matrix;

    void sanity_check()
    {
        if (mini_size >= kmer_size) {
            throw std::runtime_error("minimizer size must be smaller than k-mer size");
        }
    }
};

using muset_pa_options_t = std::shared_ptr<struct muset_pa_options>;


class musetPaCli
{

public:
    musetPaCli(
        const std::string& name,
        const std::string& desc,
        const std::string& version,
        const std::string& authors
    );

    muset_pa_options_t parse(int argc, char* argv[]);

private:
    muset_pa_cli_t cli{nullptr};
    muset_pa_options_t muset_pa_opt{nullptr};
};



muset_pa_options_t muset_pa_cli(muset_pa_cli_t cli, muset_pa_options_t options);

};  // namespace km
