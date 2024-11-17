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
  filter_opt = std::make_shared<struct filter_options>();
  merge_opt = std::make_shared<struct merge_options>();

  fafmt_cli(cli, fafmt_opt);
  fasta_cli(cli, fasta_opt);
  filter_cli(cli, filter_opt);
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
    }
    else if (cli->is("fasta")) {
        this->fasta_opt->inputs = cli->get_positionals();
        return std::make_tuple(COMMAND::FASTA, this->fasta_opt);
    }
    else if (cli->is("merge")) {
        this->merge_opt->inputs = cli->get_positionals();
        return std::make_tuple(COMMAND::MERGE, this->merge_opt);
    }
    else if (cli->is("filter")) {
        this->filter_opt->inputs = cli->get_positionals();
        return std::make_tuple(COMMAND::FILTER, this->filter_opt);
    }
    else {
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

    fafmt->set_positionals(1, "<input.fa>", "<input.fa> - FASTA file");

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

    fasta->set_positionals(1, "<in.mat>", "<in.mat> - text-based k-mer matrix");

    return opt;
}


kmat_opt_t filter_cli(std::shared_ptr<bc::Parser<1>> cli, filter_opt_t opt)
{
    bc::cmd_t filter = cli->add_command("filter", "Filter a matrix by selecting k-mers that are potentially differential.");

    filter->add_group("main options", "");

    filter->add_param("-o/--output", "output file.")
        ->meta("FILE")
        ->def("")
        ->setter(opt->output);

    // filter->add_param("--cpr", "enable compressed binary output.")
    //     ->as_flag()
    //     ->setter(opt->lz4);

    filter->add_group("filtering options", "");

    filter->add_param("-a/--min-abundance", "min abundance to keep a k-mer.")
        ->meta("INT")
        ->def("1")
        ->checker(bc::check::is_number)
        ->setter(opt->min_abundance);

    filter->add_param("-f/--min-frac-absent", "fraction of samples from which a k-mer should be absent. [0.0, 1.0]")
        ->meta("FLOAT")
        ->def("0.1")
        ->checker(bc::check::f::range(0.0, 1.0))
        ->setter(opt->min_frac_absent);

    filter->add_param("-F/--min-frac-present", "fraction of samples in which a k-mer should be present. [0.0, 1.0]")
        ->meta("FLOAT")
        ->def("0.1")
        ->checker(bc::check::f::range(0.0, 1.0))
        ->setter(opt->min_frac_present);

    filter->add_param("-n/--min-nb-absent", "minimum number of samples from which a k-mer should be absent (overrides -f).")
        ->meta("FLOAT")
        ->def("0")
        ->checker(bc::check::is_number)
        ->setter(opt->min_nb_absent)
        ->callback([opt](){ opt->min_nb_absent_set = true; });

    filter->add_param("-N/--min-nb-present", "minimum number of samples in which a k-mer should be present (overrides -F).")
        ->meta("FLOAT")
        ->def("0")
        ->checker(bc::check::is_number)
        ->setter(opt->min_nb_present)
        ->callback([opt](){ opt->min_nb_present_set = true; });

    filter->add_group("other options", "");
    
    filter->add_param("--keep-tmp", "keep temporary files.")
        ->as_flag()
        ->setter(opt->keep_tmp);

    filter->add_param("-t/--threads", "number of threads.")
        ->meta("INT")
        ->def("4")
        ->checker(bc::check::is_number)
        ->setter(opt->nb_threads);

    filter->add_param("-h/--help", "show this message and exit.")
         ->as_flag()
         ->action(bc::Action::ShowHelp);
    
    filter->add_param("-v/--version", "show version and exit.")
         ->as_flag()
         ->action(bc::Action::ShowVersion);

    filter->set_positionals(1, "<INPUT>", "<INPUT> - text-based k-mer matrix or kmtricks run directory");

    return opt;
}


kmat_opt_t merge_cli(std::shared_ptr<bc::Parser<1>> cli, merge_opt_t opt)
{
    bc::cmd_t merge = cli->add_command("merge", "Merge two input text-based kmer-sorted matrices.");

    merge->add_param("-k/--kmer-size", "k-mer size")
         ->meta("INT")
         ->def("31")
         ->setter(opt->kmer_size);
    
    merge->add_param("-o/--output", "output file.")
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

    merge->set_positionals(2, "<MATRIX_1> <MATRIX_2>", "Two text-based k-mer matrices");

    return opt;
}

};