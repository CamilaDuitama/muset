#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <kmat_tools/cmd/fasta.h>
#include <kmat_tools/matrix.h>
#include <kmat_tools/utils.h>


namespace kmat {

int main_fasta(fasta_opt_t opt)
{
    // opt->sanity_check();

    TextMatrixReader mat(opt->inputs[0]);

    std::ostream* fpout = &std::cout;
    std::ofstream ofs;
    if(!(opt->output).empty()) {
        ofs.open((opt->output).c_str());
        if(!ofs.good()) { throw std::runtime_error(fmt::format("cannot open {}", opt->output)); }
        fpout = &ofs;
    }

    std::string kmer;
    size_t kmer_count{0};
    while (mat.read_kmer(kmer)) {
        bool valid_kmer = kmer.length() > 0 && std::all_of(kmer.begin(), kmer.end(), [](const char c) { return isnuc[c]; });
        if(!valid_kmer) {
            spdlog::warn(fmt::format("skipping invalid k-mer at line {}: \"{}\"", mat.line_count(), kmer));
            continue;
        }
        *fpout << ">" << ++kmer_count << "\n";
        *fpout << kmer << "\n";
    }

    if(!(opt->output).empty()) {
        ofs.close();
    }

    spdlog::info(fmt::format("{} k-mers processed", kmer_count));

    return 0;
}

};