#pragma once

#include <memory>
#include <string>
#include <vector>

#include <fmt/format.h>

#include <kmat_tools/cli/cli_common.h>
#include <kmat_tools/cli/reverse.h>
#include <kmat_tools/matrix.h>
#include <kmat_tools/utils.h>

namespace kmat {

int main_reverse(reverse_opt_t opt)
{
    // opt->sanity_check();

    TextMatrixReader reader(opt->inputs[0]);

    std::ostream* fpout = &std::cout;
    std::ofstream ofs;
    if(!(opt->output).empty()) {
        ofs.open((opt->output).c_str());
        if(!ofs.good()) { throw std::runtime_error(fmt::format("cannot open {}", opt->output)); }
        fpout = &ofs;
    }

    std::string kmer, line;
    while (reader.read_kmer_and_line(kmer,line)) {
        if(!opt->canonicalize || !is_canonical(kmer,opt->actg_order)) {
            kmat::reverse_complement_inplace(kmer);
        }
        *fpout << kmer << " " << line << "\n";
    }

    if(!(opt->output).empty()) {
        ofs.close();
    }

    return 0;
}

};