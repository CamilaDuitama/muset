#include <filesystem>
#include <kmtricks/utils.hpp>
#include <fmt/format.h>

#include "muset_pa_cli.h"

namespace fs = std::filesystem;


namespace muset {

musetPaCli::musetPaCli(const std::string& name, const std::string& desc, const std::string& version, const std::string& authors) {
    cli = std::make_shared<bc::Parser<0>>(name, desc, version, authors);
    muset_pa_opt = std::make_shared<struct muset_pa_options>();
    muset_pa_cli(cli, muset_pa_opt);
}

muset_pa_options_t musetPaCli::parse(int argc, char* argv[])
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

  return muset_pa_opt;
}


auto file_exists = [](const std::string& p, const std::string& v) {
    return std::make_tuple(fs::is_regular_file(v),
        bc::utils::format_error(p, v, "File does not exist!"));
};

auto dir_already_exists = [](const std::string& p, const std::string& v) {
    return std::make_tuple(!fs::is_directory(v),
        bc::utils::format_error(p, v, "Directory already exists!"));
};


muset_pa_options_t muset_pa_cli(muset_pa_cli_t cli, muset_pa_options_t options) {

    cli->add_group("main options", "");

    cli->add_param("--file", "kmtricks-like input file, see README.md.")
       ->meta("FILE")
       ->def("")
       ->checker(file_exists)
       ->setter(options->fof);

    cli->add_param("-k/--kmer-size", "k-mer size. [8, 63].")
        ->meta("INT")
        ->def("31")
        ->checker(bc::check::f::range(8, 63))
        ->setter(options->kmer_size);

    cli->add_param("-o/--output", "output directory.")
        ->meta("DIR")
        ->def("output_pa")
        ->checker(dir_already_exists)
        ->setter(options->out_dir);

    cli->add_param("-a/--min-abundance", "minimum abundance required to keep a kmer.")
        ->meta("INT")
        ->def("2")
        ->checker(bc::check::is_number)
        ->setter(options->min_abundance);

    cli->add_param("-l/--min-unitig-length", "minimum length required to keep a unitig.")
        ->meta("INT")
        ->def("2k-1")
        ->setter(options->min_utg_len)
        ->callback([options](){ options->min_utg_len_set = true; });

    cli->add_param("-r/--min-utg-frac", "output a binary matrix, using this minimum treshold to set a unitig as present (1) in a sample [0,1].")
        ->meta("FLOAT")
        ->def("0.8")
        ->checker(bc::check::f::range(0.0, 1.0))
        ->setter(options->min_utg_frac)
        ->callback([options](){ options->min_utg_frac_set = true; });

    cli->add_param("-s/--write-seq", "write the unitig sequence instead of the identifier in the output matrix")
        ->as_flag()
        ->setter(options->write_utg_seq);

    cli->add_param("-m/--mini-size", "minimizer size. [4, 15].")
        ->meta("INT")
        ->def("15")
        ->checker(bc::check::f::range(4, 15))
        ->setter(options->mini_size);

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
