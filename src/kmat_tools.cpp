#include <cstdlib>
#include <memory>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <kmat_tools/cmd.h>
#include <kmat_tools/cli.h>


int main(int argc, char* argv[])
{
    kmat::kmatCli cli("kmat_tools", "a collection of tools to process text-based k-mer matrices", PROJECT_VER, "");
    auto [cmd, options] = cli.parse(argc, argv);

    // set_verbosity_level(options->verbosity);
    auto cerr_logger = spdlog::stderr_color_mt("kmat_tools");
    cerr_logger->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] %v");
    spdlog::set_default_logger(cerr_logger);

    spdlog::info("Command: {}", cmd_to_str(cmd));

    try
    {
        if (cmd == kmat::COMMAND::CONVERT) {
            kmat::convert_opt_t opt = std::static_pointer_cast<struct kmat::convert_options>(options);
            return kmat::main_convert(opt);
        }
        else if (cmd == kmat::COMMAND::DIFF) {
            kmat::diff_opt_t opt = std::static_pointer_cast<struct kmat::diff_options>(options);
            return kmat::main_diff(opt);
        }
        else if (cmd == kmat::COMMAND::FAFMT) {
            kmat::fafmt_opt_t opt = std::static_pointer_cast<struct kmat::fafmt_options>(options);
            return kmat::main_fafmt(opt);
        }
        else if (cmd == kmat::COMMAND::FAFMT) {
            kmat::fafmt_opt_t opt = std::static_pointer_cast<struct kmat::fafmt_options>(options);
            return kmat::main_fafmt(opt);
        }
        else if (cmd == kmat::COMMAND::FASTA) {
            kmat::fasta_opt_t opt = std::static_pointer_cast<struct kmat::fasta_options>(options);
            return kmat::main_fasta(opt);
        }
        else if (cmd == kmat::COMMAND::FILTER) {
            kmat::filter_opt_t opt = std::static_pointer_cast<struct kmat::filter_options>(options);
            return kmat::main_filter(opt);
        }
        else if (cmd == kmat::COMMAND::MERGE) {
            kmat::merge_opt_t opt = std::static_pointer_cast<struct kmat::merge_options>(options);
            return kmat::main_merge(opt);
        }
        else if (cmd == kmat::COMMAND::REVERSE) {
            kmat::reverse_opt_t opt = std::static_pointer_cast<struct kmat::reverse_options>(options);
            return kmat::main_reverse(opt);
        }
    }
    catch (const std::exception& e)
    {
        spdlog::error(e.what());
        exit(EXIT_FAILURE);
    }

    return EXIT_SUCCESS;
}