#pragma once

#include <memory>
#include <string>

#include <algorithm>
#include <bcli/bcli.hpp>

#include <kmat_tools/config.h>
#include <kmat_tools/cli/cli_common.h>
#include <kmat_tools/cli/convert.h>
#include <kmat_tools/cli/diff.h>
#include <kmat_tools/cli/fafmt.h>
#include <kmat_tools/cli/fasta.h>
#include <kmat_tools/cli/filter.h>
#include <kmat_tools/cli/merge.h>

namespace kmat
{

class kmatCli
{

public:
  kmatCli(
      const std::string& name,
      const std::string& desc,
      const std::string& version,
      const std::string& authors);

  std::tuple<COMMAND, kmat_opt_t> parse(int argc, char* argv[]);

private:
  
  cli_t cli {nullptr};
  
  convert_opt_t convert_opt {nullptr};
  diff_opt_t    diff_opt {nullptr};
  fafmt_opt_t   fafmt_opt {nullptr};
  fasta_opt_t   fasta_opt {nullptr};
  filter_opt_t  filter_opt {nullptr};
  merge_opt_t   merge_opt {nullptr};  
};

};  // namespace kmat