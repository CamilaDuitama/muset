#include <memory>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <kmat_tools/cmd/select.h>
#include <kmat_tools/matrix.h>
#include <kmat_tools/utils.h>


namespace kmat {

int main_select(select_opt_t opt)
{
    TextMatrixReader mat_1(opt->inputs[0]);
    TextMatrixReader mat_2(opt->inputs[1]);

    std::ostream* fpout = &std::cout;
    std::ofstream ofs;
    if(!(opt->output).empty()) {
        ofs.open((opt->output).c_str());
        if(!ofs.good()) { throw std::runtime_error(fmt::format("cannot open {}", opt->output)); }
        fpout = &ofs;
    }

    std::string kmer_1;
    bool has_kmer_1 = mat_1.read_kmer(kmer_1);

    std::string kmer_2, line_2;
    bool has_kmer_2 = mat_2.read_kmer_and_line(kmer_2, line_2);

    if(has_kmer_1 && has_kmer_2 && kmer_1.size() != kmer_2.size()) {
        throw std::runtime_error(fmt::format("different k-mer size between the two input matrices: {} vs {}.", kmer_1.size(), kmer_2.size()));
    }

    size_t retained_kmers{0};
    while (has_kmer_1 && has_kmer_2) {
        int cmp = opt->actg_order ? actg_compare(kmer_1,kmer_2) : kmer_1.compare(kmer_2);
        if(cmp == 0) { // kmer_1 == kmer_2
            retained_kmers++;
            *fpout << kmer_2 << " " << line_2 << "\n";
            has_kmer_1 = mat_1.read_kmer(kmer_1);
            has_kmer_2 = mat_2.read_kmer_and_line(kmer_2, line_2);
        } else if(cmp < 0) { // kmer_1 < kmer_2
            has_kmer_1 = mat_1.read_kmer(kmer_1);
        } else { // kmer_1 > kmer_2
            has_kmer_2 = mat_2.read_kmer_and_line(kmer_2, line_2);
        }
    }

    spdlog::info(fmt::format("retained k-mers: {}", retained_kmers));

    if(!(opt->output).empty()) {
        ofs.close();
    }

    return 0;
}

};