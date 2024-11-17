#pragma once

#include <kmat_tools/cli/cli_common.h>

namespace kmat {

struct fasta_options : kmat_options
{
  std::string input;
  std::string output;
};

using fasta_opt_t = std::shared_ptr<struct fasta_options>;

kmat_opt_t fasta_cli(std::shared_ptr<bc::Parser<1>> cli, fasta_opt_t options);

};