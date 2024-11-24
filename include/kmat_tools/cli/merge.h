#pragma once

#include <kmat_tools/cli/cli_common.h>

namespace kmat {

struct merge_options : kmat_options
{
  uint32_t kmer_size{31};
  std::string output;
  bool actg_order{false};
};

using merge_opt_t = std::shared_ptr<struct merge_options>;

kmat_opt_t merge_cli(std::shared_ptr<bc::Parser<1>> cli, merge_opt_t options);

};