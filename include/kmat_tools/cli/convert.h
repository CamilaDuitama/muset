#pragma once

#include <kmat_tools/cli/cli_common.h>

namespace kmat {

struct convert_options : kmat_options
{
    std::string out_fname;

    bool ap_flag{false};
    double min_frac{0.8};
    bool min_frac_set{false};
    bool out_write_seq{false};
    bool no_header{false};
    bool out_csv{false};
};

using convert_opt_t = std::shared_ptr<struct convert_options>;

kmat_opt_t convert_cli(std::shared_ptr<bc::Parser<1>> cli, convert_opt_t options);

};