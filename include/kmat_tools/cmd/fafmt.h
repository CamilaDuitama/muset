#pragma once

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <kseq++/seqio.hpp>

#include <kmat_tools/cli/cli_common.h>
#include <kmat_tools/cli/fafmt.h>
#include <kmat_tools/utils.h>

namespace kmat {

int main_fafmt(fafmt_opt_t opt)
{
    std::string input_file = opt->inputs[0];
    if(!std::filesystem::exists(input_file)) {
        throw std::runtime_error(fmt::format("input file {} does not exist", input_file));
    }

    std::ostream* fpout = &std::cout;
    std::ofstream ofs;
    if(!(opt->output).empty()) {
        ofs.open((opt->output).c_str());
        if(!ofs.good()) { throw std::runtime_error(fmt::format("cannot open output file {}", opt->output)); }
        fpout = &ofs;
    }

    klibpp::KSeq record;
    klibpp::SeqStreamIn ssi(input_file.c_str());
    size_t total{0};
    size_t retained{0};
    while (ssi >> record) {
        total++;
        if (record.seq.length() < opt->min_length) { continue; }
        *fpout << '>' << record.name;
        if(!record.comment.empty()) { *fpout << ' ' << record.comment; }
        *fpout << '\n' << record.seq << '\n';
        retained++;
    }

    if(!(opt->output).empty()) {
        ofs.close();
    }

    spdlog::info(fmt::format("{}/{} sequences retained", retained, total));

    return 0;
}

};