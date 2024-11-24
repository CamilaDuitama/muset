#pragma once

#include <filesystem>

#define WITH_KM_IO
#include <kmtricks/public.hpp>

#include <kmat_tools/cli/cli_common.h>

namespace fs = std::filesystem;

namespace kmat {

struct filter_options : kmat_options {

    fs::path output;
    
    uint32_t min_abundance{0};

    double min_frac_absent{0.1};
    double min_frac_present{0.1};
    bool min_nb_absent_set{false};
    int min_nb_absent{0};
    bool min_nb_present_set{false};
    int min_nb_present{0};
    
    bool keep_tmp{false};
    bool lz4{true}; // not used yet

    // kmtricks-related options
    fs::path matrices_dir;
    fs::path filtered_dir;
    uint32_t kmer_size{31};

    size_t nb_threads{1};
};

using filter_opt_t = std::shared_ptr<struct filter_options>;

kmat_opt_t filter_cli(std::shared_ptr<bc::Parser<1>> cli, filter_opt_t options);

};