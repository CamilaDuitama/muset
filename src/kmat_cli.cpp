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
  
    convert_opt = std::make_shared<struct convert_options>();
    diff_opt = std::make_shared<struct diff_options>();
    fafmt_opt = std::make_shared<struct fafmt_options>();
    fasta_opt = std::make_shared<struct fasta_options>();
    filter_opt = std::make_shared<struct filter_options>();
    merge_opt = std::make_shared<struct merge_options>();
    reverse_opt = std::make_shared<struct reverse_options>();
    select_opt = std::make_shared<struct select_options>();
    unitig_opt = std::make_shared<struct unitig_options>();

    convert_cli(cli, convert_opt);
    diff_cli(cli, diff_opt);
    fafmt_cli(cli, fafmt_opt);
    fasta_cli(cli, fasta_opt);
    filter_cli(cli, filter_opt);
    merge_cli(cli, merge_opt);
    reverse_cli(cli, reverse_opt);
    unitig_cli(cli, unitig_opt);
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
        fmt::print("kmtricks {}\n", MUSET_PROJECT_VER);
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

    if (cli->is("convert")) {
        this->convert_opt->inputs = cli->get_positionals();
        return std::make_tuple(COMMAND::CONVERT, this->convert_opt);
    }
    else if (cli->is("diff")) {
        this->diff_opt->inputs = cli->get_positionals();
        return std::make_tuple(COMMAND::DIFF, this->diff_opt);
    }
    else if (cli->is("fafmt")) {
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
    else if (cli->is("reverse")) {
        this->reverse_opt->inputs = cli->get_positionals();
        return std::make_tuple(COMMAND::REVERSE, this->reverse_opt);
    }
    else if (cli->is("select")) {
        this->select_opt->inputs = cli->get_positionals();
        return std::make_tuple(COMMAND::SELECT, this->select_opt);
    }
    else if (cli->is("unitig")) {
        this->unitig_opt->inputs = cli->get_positionals();
        return std::make_tuple(COMMAND::UNITIG, this->unitig_opt);
    }
    else {
        return std::make_tuple(COMMAND::UNKNOWN, std::make_shared<struct kmat_options>());
    }
}


kmat_opt_t convert_cli(std::shared_ptr<bc::Parser<1>> cli, convert_opt_t opt)
{
    bc::cmd_t convert = cli->add_command("convert", "Convert ggcat jsonl color output into a unitig matrix");

    convert->add_param("-p/--pa-out", "output a binary presence-absence matrix, instead of k-mer fractions")
        ->as_flag()
        ->setter(opt->ap_flag);

    convert->add_param("-f/--min-frac", "minimum fraction to set presence to 1 (only effective with -p). [0.5, 1.0]")
        ->meta("FLOAT")
        ->def("0.8")
        ->checker(bc::check::f::range(0.5, 1.0))
        ->setter(opt->min_frac)
        ->callback([opt](){ opt->min_frac_set = true; });

    convert->add_param("-s/--write-seq", "write unitig sequence instead of ID in the output matrix")
        ->as_flag()
        ->setter(opt->out_write_seq);

    convert->add_param("-H/--no-header", "do not write header in output matrix")
        ->as_flag()
        ->setter(opt->no_header);

    convert->add_param("--csv", "write csv output")
        ->as_flag()
        ->setter(opt->out_csv);

    convert->add_param("-o/--output", "output file. {stdout}")
         ->meta("FILE")
         ->def("")
         ->setter(opt->out_fname);
    
    convert->add_param("-h/--help", "show this message and exit.")
         ->as_flag()
         ->action(bc::Action::ShowHelp);
    
    convert->add_param("-v/--version", "show version and exit.")
         ->as_flag()
         ->action(bc::Action::ShowVersion);

    convert->set_positionals(3,
        "<unitigs.fasta> <color_names_dump.jsonl> <query_output.jsonl>",
        "");

    return opt;
}


kmat_opt_t diff_cli(std::shared_ptr<bc::Parser<1>> cli, diff_opt_t opt)
{
    bc::cmd_t diff = cli->add_command("diff", "Difference between two sorted k-mer matrices");
    
    diff->add_param("-o/--output", "output file. {stdout}")
         ->meta("FILE")
         ->def("")
         ->setter(opt->output);
    
    diff->add_param("-z/--actg", "use A<C<T<G order of nucleotides")
         ->as_flag()
         ->setter(opt->actg_order);

    diff->add_param("-h/--help", "show this message and exit.")
         ->as_flag()
         ->action(bc::Action::ShowHelp);
    
    diff->add_param("-v/--version", "show version and exit.")
         ->as_flag()
         ->action(bc::Action::ShowVersion);

    diff->set_positionals(2, "<matrix_1> <matrix_2>", "Two text-based k-mer matrices");

    return opt;
}


kmat_opt_t fafmt_cli(std::shared_ptr<bc::Parser<1>> cli, fafmt_opt_t opt)
{
    bc::cmd_t fafmt = cli->add_command("fafmt", "Filter a FASTA file by length and write sequences in single lines");

    fafmt->add_param("-l/--min-length", "minimum sequence length")
         ->meta("INT")
         ->def("0")
         ->setter(opt->min_length);

    fafmt->add_param("-o/--output", "output file. {stdout}")
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

    fasta->add_param("-o/--output", "output file. {stdout}")
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

    filter->add_param("-o/--output", "output file. {stdout}")
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
    
    merge->add_param("-o/--output", "output file. {stdout}")
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

    merge->set_positionals(2, "<matrix_1> <matrix_2>", "Two text-based k-mer matrices");

    return opt;
}


kmat_opt_t reverse_cli(std::shared_ptr<bc::Parser<1>> cli, reverse_opt_t opt)
{
    bc::cmd_t reverse = cli->add_command("reverse", "Reverse-complement k-mers in a k-mer matrix file.");
    
    reverse->add_param("-o/--output", "output file. {stdout}")
         ->meta("FILE")
         ->def("")
         ->setter(opt->output);

    reverse->add_param("-c/--canonical", "output canonical sequence instead of reverse complement.")
         ->as_flag()
         ->setter(opt->canonicalize);

    reverse->add_param("-z/--actg", "use A<C<T<G order of nucleotides (only effective with -c).")
         ->as_flag()
         ->setter(opt->actg_order);

    reverse->add_param("-h/--help", "show this message and exit.")
         ->as_flag()
         ->action(bc::Action::ShowHelp);
    
    reverse->add_param("-v/--version", "show version and exit.")
         ->as_flag()
         ->action(bc::Action::ShowVersion);

    reverse->set_positionals(1, "<input.mat>", "Text-based k-mer matricx");

    return opt;
}


kmat_opt_t select_cli(std::shared_ptr<bc::Parser<1>> cli, select_opt_t opt)
{
    bc::cmd_t select = cli->add_command("select", "Select from an input matrix only k-mers that belong to a reference matrix.");
    
    select->add_param("-o/--output", "output file. {stdout}")
         ->meta("FILE")
         ->def("")
         ->setter(opt->output);

    select->add_param("-z/--actg", "use A<C<T<G order of nucleotides")
         ->as_flag()
         ->setter(opt->actg_order);

    select->add_param("-h/--help", "show this message and exit.")
         ->as_flag()
         ->action(bc::Action::ShowHelp);
    
    select->add_param("-v/--version", "show version and exit.")
         ->as_flag()
         ->action(bc::Action::ShowVersion);

    select->set_positionals(2, "<reference.mat> <input.mat>", "Reference and input (sorted) text-based k-mer matrices");

    return opt;
}


kmat_opt_t unitig_cli(std::shared_ptr<bc::Parser<1>> cli, unitig_opt_t opt)
{
    bc::cmd_t unitig = cli->add_command("unitig", "Create a unitig matrix.");

    unitig->add_group("main options", "");

    unitig->add_param("-k/--kmer-size", fmt::format("k-mer size [8,{}].", KL[MUSET_KMER_N-1]-1))
        ->meta("INT")
        ->def("31")
        ->checker(bc::check::f::range(8, KL[MUSET_KMER_N-1]-1))
        ->setter(opt->kmer_size);

    unitig->add_param("-p/--prefix", "output files prefix.")
        ->meta("FILE")
        ->def("out")
        ->setter(opt->prefix);

    unitig->add_param("-f/--min-frac", "set average abundance to 0 if k-mer fraction is below this threshold [0,1].")
        ->meta("FLOAT")
        ->def("0.0")
        ->checker(bc::check::f::range(0.0, 1.0))
        ->setter(opt->min_frac);

    unitig->add_param("--out-frac", "output an additional matrix containing k-mer fractions.")
        ->as_flag()
        ->setter(opt->write_frac_matrix);

    unitig->add_param("-s/--write-seq", "write the unitig sequence instead of the identifier in the output matrix")
        ->as_flag()
        ->setter(opt->write_seq);

    unitig->add_group("other options", "");

    unitig->add_param("-m/--minimizer-size", "minimizer size")
        ->meta("INT")
        ->def("15")
        ->setter(opt->mini_size);

    unitig->add_param("-t/--threads", "number of threads.")
        ->meta("INT")
        ->def("4")
        ->checker(bc::check::is_number)
        ->setter(opt->nb_threads);

    unitig->add_param("-h/--help", "show this message and exit.")
         ->as_flag()
         ->action(bc::Action::ShowHelp);
    
    unitig->add_param("-v/--version", "show version and exit.")
         ->as_flag()
         ->action(bc::Action::ShowVersion);

    unitig->set_positionals(2, "<unitigs.fasta> <kmer_matrix>", "a unitig fasta file and a text-based k-mer matrix");

    return opt;
}

};