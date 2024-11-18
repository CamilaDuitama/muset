#pragma once

#include <kmat_tools/cli/cli_common.h>

namespace kmat {

struct diff_options : kmat_options
{
    std::string output;

    uint32_t kmer_size{31};
    bool actg_order{false};
};

using diff_opt_t = std::shared_ptr<struct diff_options>;

kmat_opt_t diff_cli(std::shared_ptr<bc::Parser<1>> cli, diff_opt_t options);

};