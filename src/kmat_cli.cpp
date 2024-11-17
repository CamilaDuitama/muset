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
  
  fafmt_opt = std::make_shared<struct fafmt_options>();
  fasta_opt = std::make_shared<struct fasta_options>();
  merge_opt = std::make_shared<struct merge_options>();

  fafmt_cli(cli, fafmt_opt);
  fasta_cli(cli, fasta_opt);
  merge_cli(cli, merge_opt);
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

    if (cli->is("fafmt")) {
        this->fafmt_opt->inputs = cli->get_positionals();
        return std::make_tuple(COMMAND::FAFMT, this->fafmt_opt);
    } else if (cli->is("fasta")) {
        this->fasta_opt->inputs = cli->get_positionals();
        return std::make_tuple(COMMAND::FASTA, this->fasta_opt);
    } else if (cli->is("merge")) {
        this->merge_opt->inputs = cli->get_positionals();
        return std::make_tuple(COMMAND::MERGE, this->merge_opt);
    } else {
        return std::make_tuple(COMMAND::UNKNOWN, std::make_shared<struct kmat_options>());
    }
}


kmat_opt_t fafmt_cli(std::shared_ptr<bc::Parser<1>> cli, fafmt_opt_t opt)
{
    bc::cmd_t fafmt = cli->add_command("fafmt", "Filter a FASTA file by length and write sequences in single lines");

    fafmt->add_param("-l/--min-length", "minimum sequence length")
         ->meta("INT")
         ->def("0")
         ->setter(opt->min_length);

    fafmt->add_param("-o/--output", "output file.")
         ->meta("FILE")
         ->def("")
         ->setter(opt->output);
    
    fafmt->add_param("-h/--help", "show this message and exit.")
         ->as_flag()
         ->action(bc::Action::ShowHelp);
    
    fafmt->add_param("-v/--version", "show version and exit.")
         ->as_flag()
         ->action(bc::Action::ShowVersion);

    fafmt->set_positionals(1, "<input.fa>", "");

    return opt;
}


kmat_opt_t fasta_cli(std::shared_ptr<bc::Parser<1>> cli, fasta_opt_t opt)
{
    bc::cmd_t fasta = cli->add_command("fasta", "Output a k-mer matrix in FASTA format");

    fasta->add_param("-o/--output", "output file.")
         ->meta("FILE")
         ->def("")
         ->setter(opt->output);
    
    fasta->add_param("-h/--help", "show this message and exit.")
         ->as_flag()
         ->action(bc::Action::ShowHelp);
    
    fasta->add_param("-v/--version", "show version and exit.")
         ->as_flag()
         ->action(bc::Action::ShowVersion);

    fasta->set_positionals(1, "<in.mat>", "");

    return opt;
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