#pragma once

#include <memory>
#include <string>
#include <vector>

#include <fmt/format.h>

#include <kmat_tools/cli/cli_common.h>
#include <kmat_tools/cli/diff.h>
#include <kmat_tools/matrix.h>
#include <kmat_tools/utils.h>

namespace kmat {

int main_diff(diff_opt_t opt)
{
  // opt->sanity_check();

  TextMatrixReader mat_1(opt->inputs[0]);
  TextMatrixReader mat_2(opt->inputs[1]);

  std::ostream* fpout = &std::cout;
  std::ofstream ofs;
  if(!(opt->output).empty()) {
    ofs.open((opt->output).c_str());
    if(!ofs.good()) { throw std::runtime_error(fmt::format("cannot open {}", opt->output)); }
    fpout = &ofs;
  }

  std::string kmer_1, line_1;
  bool has_kmer_1 = mat_1.read_kmer_and_line(kmer_1, line_1);

  std::string kmer_2, line_2;
  bool has_kmer_2 = mat_2.read_kmer_and_line(kmer_2, line_2);

  if(has_kmer_1 && has_kmer_2 && kmer_1.size() != kmer_2.size()) {
    throw std::runtime_error(fmt::format("different k-mer size between the two input matrices: {} vs {}.", kmer_1.size(), kmer_2.size()));
  }

  while (has_kmer_1 && has_kmer_2) {
    int cmp = opt->actg_order ? actg_compare(kmer_1,kmer_2) : kmer_1.compare(kmer_2);
    if(cmp == 0) { // kmer_1 == kmer_2
      has_kmer_1 = mat_1.read_kmer_and_line(kmer_1, line_1);
      has_kmer_2 = mat_2.read_kmer_and_line(kmer_2, line_2);
    } else if(cmp < 0) { // kmer_1 < kmer_2
      *fpout << kmer_1 << " " << line_1 << "\n";
      has_kmer_1 = mat_1.read_kmer_and_line(kmer_1, line_1);
    } else { // kmer_1 > kmer_2
      has_kmer_2 = mat_2.read_kmer_and_line(kmer_2, line_2);
    }
  }

  while (has_kmer_1) {
    *fpout << kmer_1 << " " << line_1 << "\n";
    has_kmer_1 = mat_1.read_kmer_and_line(kmer_1, line_1);
  }

  if(!(opt->output).empty()) {
    ofs.close();
  }

  return 0;
}

};