#pragma once

#include <kmat_tools/cli/cli_common.h>

namespace kmat {

struct select_options : kmat_options
{
    std::string output;

    uint32_t kmer_size{31};
    bool actg_order{false};
};

using select_opt_t = std::shared_ptr<struct select_options>;

kmat_opt_t select_cli(std::shared_ptr<bc::Parser<1>> cli, select_opt_t options);

};