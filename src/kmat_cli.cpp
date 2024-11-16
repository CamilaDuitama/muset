#include <kmat_tools/cli.h>
#include <filesystem>

#include <fmt/format.h>

namespace fs = std::filesystem;
namespace kmat {

kmatCli::kmatCli(
    const std::string& name,
    const std::string& desc,
    const std::string& version,
    const std::string& authors)
{
  cli = std::make_shared<bc::Parser<1>>(bc::Parser<1>(name, desc, version, authors));
  
  merge_opt = std::make_shared<struct merge_options>();
  //filter_opt = std::make_shared<struct filter_options>(filter_options{});

  merge_cli(cli, merge_opt);
  //filter_cli(cli, filter_opt);
}

std::tuple<COMMAND, kmat_opt_t> kmatCli::parse(int argc, char* argv[])
{
    if (argc > 1 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h"))
    {
        cli->show_help();
        exit(EXIT_FAILURE);
    }
    else if (argc > 1 && (std::string(argv[1]) == "--version" || std::string(argv[1]) == "-v"))
    {
        fmt::print("kmtricks {}\n", PROJECT_VER);
        exit(EXIT_FAILURE);
    }

    try
    {
        (*cli).parse(argc, argv);
    }
    catch (const bc::ex::BCliError& e)
    {
        bc::utils::exit_bcli(e);
        exit(EXIT_FAILURE);
    }

    if (cli->is("merge")) {
        this->merge_opt->inputs = cli->get_positionals();
        return std::make_tuple(COMMAND::MERGE, this->merge_opt);
    } else {
        return std::make_tuple(COMMAND::UNKNOWN, std::make_shared<struct kmat_options>());
    }
}


kmat_opt_t merge_cli(std::shared_ptr<bc::Parser<1>> cli, merge_opt_t opt)
{
    bc::cmd_t merge = cli->add_command("merge", "Merge two input text-based kmer-sorted matrices.");

    merge->add_param("-k/--kmer-size", "k-mer size")
         ->meta("INT")
         ->def("31")
         ->setter(opt->kmer_size);
    
    merge->add_param("-o/--output", "output directory.")
         ->meta("FILE")
         ->def("")
         ->setter(opt->output);
    
    merge->add_param("-z/--actg", "use A<C<T<G order of nucleotides")
         ->as_flag()
         ->setter(opt->actg_order);

    merge->add_param("-h/--help", "show this message and exit.")
         ->as_flag()
         ->action(bc::Action::ShowHelp);
    
    merge->add_param("-v/--version", "show version and exit.")
         ->as_flag()
         ->action(bc::Action::ShowVersion);

    merge->set_positionals(2, "<matrix_1> <matrix_2>", "");

    return opt;
}

};