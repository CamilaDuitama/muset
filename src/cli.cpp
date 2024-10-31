#include <filesystem>
#include <kmtricks/utils.hpp>
#include <fmt/format.h>

#include "common.h"
#include "cli.h"

namespace fs = std::filesystem;
namespace muset {

musetCli::musetCli(const std::string& name, const std::string& desc, const std::string& version, const std::string& authors) {
    cli = std::make_shared<bc::Parser<0>>(name, desc, version, authors);
    muset_opt = std::make_shared<struct muset_options>(muset_options{});
    muset_cli(cli, muset_opt);
}


muset_options_t musetCli::parse(int argc, char* argv[])
{

  try
  {
    (*cli).parse(argc, argv);
  }
  catch (const bc::ex::BCliError& e)
  {
    bc::utils::exit_bcli(e);
    exit(EXIT_FAILURE);
  }

  return muset_opt;
}


auto file_exists = [](const std::string& p, const std::string& v) {
  return std::make_tuple( fs::is_regular_file(v),
    bc::utils::format_error(p, v, "File does not exist!") );
};

auto is_km_dir = [](const std::string& p, const std::string& v) -> bc::check::checker_ret_t {
  std::string c1 = fmt::format("{}/{}", v, "kmtricks.fof");
  std::string c2 = fmt::format("{}/{}", v, "run_infos.txt");
  return std::make_tuple(
    fs::is_regular_file(c1) && fs::is_regular_file(c2),
    bc::utils::format_error(p, v, "Not a kmtricks directory!")
  );
};


muset_options_t muset_cli(std::shared_ptr<bc::Parser<0>> cli, muset_options_t options) {

    cli->add_group("main options", "");

    cli->add_param("-i/--in-matrix", "input text matrix")
       ->meta("FILE")
       ->def("")
       ->setter(options->in_matrix);

    cli->add_param("-o/--output", "output directory.")
       ->meta("DIR")
       ->def("output")
       ->setter(options->out_dir);

    cli->add_param("-k/--kmer-size", "k-mer size. [8, 63].")
       ->meta("INT")
       ->def("31")
       ->checker(bc::check::f::range(8, 63))
       ->setter(options->kmer_size);

    cli->add_param("-m/--mini-size", "minimizer size. [4, 15].")
       ->meta("INT")
       ->def("15")
       ->checker(bc::check::f::range(4, 15))
       ->setter(options->mini_size);

    cli->add_param("-a/--min-abundance", "min abundance to keep a k-mer.")
       ->meta("INT")
       ->def("2")
       ->checker(bc::check::is_number)
       ->setter(options->min_abundance);

    cli->add_param("-l/--min-unitig-length", "min length to keep a unitig in the output matrix.")
       ->meta("INT")
       ->def("2k-1")
       ->callback([options](){ options->min_utg_len_set = true; });

    cli->add_param("-u/--logan", "input samples consist of Logan unitigs (i.e., with abundance).")
       ->as_flag()
       ->setter(options->logan);

    /*** FILTERING OPTIONS ***/
    
    cli->add_group("filtering options", "");

    cli->add_param("--no-kmer-filter", "disable filtering of k-mer matrix rows before unitig construction.")
       ->as_flag()
       ->setter(options->no_kmer_filter);

    cli->add_param("-f/--min-frac-absent", "fraction of samples from which a k-mer should be absent. [0.0, 1.0]")
       ->meta("FLOAT")
       ->def("0.1")
       ->checker(bc::check::f::range(0.0, 1.0))
       ->setter(options->min_frac_absent);

    cli->add_param("-F/--min-frac-present", "fraction of samples in which a k-mer should be present. [0.0, 1.0]")
       ->meta("FLOAT")
       ->def("0.1")
       ->checker(bc::check::f::range(0.0, 1.0))
       ->setter(options->min_frac_present);

    cli->add_param("-n/--min-nb-absent", "minimum number of samples from which a k-mer should be absent (overrides -f).")
       ->meta("FLOAT")
       ->def("0")
       ->checker(bc::check::is_number)
       ->setter(options->min_nb_absent)
       ->callback([options](){ options->min_nb_absent_set = true; });

    cli->add_param("-N/--min-nb-present", "minimum number of samples in which a k-mer should be present (overrides -F).")
       ->meta("FLOAT")
       ->def("0")
       ->checker(bc::check::is_number)
       ->setter(options->min_nb_present)
       ->callback([options](){ options->min_nb_present_set = true; });

    /*** OTHER OPTIONS ***/

    cli->add_group("other options", "");

    cli->add_param("--keep-tmp", "keep temporary files.")
       ->as_flag()
       ->setter(options->keep_tmp);

    cli->add_param("-t/--threads", "number of threads.")
       ->meta("INT")
       ->def("4")
       ->checker(bc::check::is_number)
       ->setter(options->nb_threads);

    cli->add_param("-h/--help", "show this message and exit.")
       ->as_flag()
       ->action(bc::Action::ShowHelp);
    
    cli->add_param("-v/--version", "show version and exit.")
       ->as_flag()
       ->action(bc::Action::ShowVersion);
    
    // cli->add_param("--verbose", "verbosity level [debug|info|warning|error].")
    //    ->meta("STR")
    //    ->def("info")
    //    ->checker(bc::check::f::in("debug|info|warning|error"))
    //    ->setter(options->verbosity);

  return options;
}


};
