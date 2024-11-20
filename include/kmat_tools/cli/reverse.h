#pragma once

#include <kmat_tools/cli/cli_common.h>

namespace kmat {

struct reverse_options : kmat_options
{
  std::string output;

  bool actg_order{false};
  bool canonicalize{false};
};

using reverse_opt_t = std::shared_ptr<struct reverse_options>;

kmat_opt_t reverse_cli(std::shared_ptr<bc::Parser<1>> cli, reverse_opt_t options);

};