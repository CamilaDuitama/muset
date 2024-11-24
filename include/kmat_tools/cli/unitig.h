#pragma once

#include <filesystem>

#define WITH_KM_IO
#include <kmtricks/public.hpp>

#include <kmat_tools/cli/cli_common.h>

namespace fs = std::filesystem;

namespace kmat {

struct unitig_options : kmat_options {

    uint32_t kmer_size{31};
    uint32_t mini_size{15};
    
    std::string prefix;

    double min_frac{0.0};
    bool write_seq{false};
    bool write_frac_matrix{false};
    size_t nb_threads{1};
};

using unitig_opt_t = std::shared_ptr<struct unitig_options>;

kmat_opt_t unitig_cli(std::shared_ptr<bc::Parser<1>> cli, unitig_opt_t options);

};