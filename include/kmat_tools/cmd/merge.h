#pragma once

#include <memory>
#include <string>
#include <vector>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <kmat_tools/cli/cli_common.h>
#include <kmat_tools/cli/merge.h>
#include <kmat_tools/matrix.h>
#include <kmat_tools/utils.h>

namespace kmat {

int main_merge(merge_opt_t opt)
{
    // opt->sanity_check();

    TextMatrixReader mat_1(opt->inputs[0]);
    TextMatrixReader mat_2(opt->inputs[1]);

    std::ostream* fpout = &std::cout;
    std::ofstream ofs;
    if(!(opt->output).empty()) {
        ofs.open((opt->output).c_str());
        if(!ofs.good()) { throw std::runtime_error(fmt::format("cannot to open {}", opt->output)); }
        fpout = &ofs;
    }

    size_t ksize = opt->kmer_size;

    std::string kmer_1, line_1;
    bool has_kmer_1 = mat_1.read_kmer_and_line(kmer_1, line_1);
    size_t nb_samples_1 = has_kmer_1 ? get_nb_samples(line_1) : 0;
    spdlog::info(fmt::format("samples in 1st matrix: {}", nb_samples_1));

    std::string empty_samples_1; empty_samples_1.reserve(2*nb_samples_1);
    for(auto i=0; i<nb_samples_1; ++i) { empty_samples_1.append(" 0"); }

    std::string kmer_2, line_2;
    bool has_kmer_2 = mat_2.read_kmer_and_line(kmer_2, line_2);
    size_t nb_samples_2 = has_kmer_2 ? get_nb_samples(line_2) : 0;
    spdlog::info(fmt::format("samples in 2nd matrix: {}", nb_samples_2));

    std::string empty_samples_2; empty_samples_2.reserve(2*nb_samples_2);
    for(auto i=0; i<nb_samples_2; ++i) { empty_samples_2.append(" 0"); }

    while (has_kmer_1 && has_kmer_2) {
        int cmp = opt->actg_order ? actg_compare(kmer_1,kmer_2) : kmer_1.compare(kmer_2);
        if(cmp == 0) {
            *fpout << kmer_1 << " " << line_1 << " " << line_2;
            has_kmer_1 = mat_1.read_kmer_and_line(kmer_1, line_1);
            has_kmer_2 = mat_2.read_kmer_and_line(kmer_2, line_2);
        } else if(cmp < 0) {
            *fpout << kmer_1 << " " << line_1 << empty_samples_2;
            has_kmer_1 = mat_1.read_kmer_and_line(kmer_1, line_1);
        } else { // cmp > 0
            *fpout << kmer_2 << empty_samples_1 << " " << line_2;
            has_kmer_2 = mat_2.read_kmer_and_line(kmer_2, line_2);
        }
        *fpout << "\n";
    }

    while (has_kmer_1) {
        *fpout << kmer_1 << " " << line_1 << empty_samples_2;
        has_kmer_1 = mat_1.read_kmer_and_line(kmer_1, line_1);
    }

    while(has_kmer_2) {
        *fpout << kmer_2 << empty_samples_1 << " " << line_2;
            has_kmer_2 = mat_2.read_kmer_and_line(kmer_2, line_2);
    }

    if(!(opt->output).empty()) {
        ofs.close();
    }

    return 0;
}

};