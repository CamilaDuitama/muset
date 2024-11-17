#pragma once

#include <kmat_tools/cli/cli_common.h>

namespace kmat {

struct fafmt_options : kmat_options
{
  std::string output;
  size_t min_length;
};

using fafmt_opt_t = std::shared_ptr<struct fafmt_options>;

kmat_opt_t fafmt_cli(std::shared_ptr<bc::Parser<1>> cli, fafmt_opt_t options);

};