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
    if (cmd == kmat::COMMAND::MERGE) {
      kmat::merge_opt_t opt = std::static_pointer_cast<struct kmat::merge_options>(options);
      return kmat::main_merge(opt);
    }
  }
  catch (const std::exception& e)
  {
    spdlog::error(e.what());
    exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}