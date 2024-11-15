#include<bcli/bcli.hpp>

using namespace bc;

int main(int argc, char* argv[])
{
  conf::get().help(true)
             .verbose(false)
             .debug(false)
             .version(true)
             .default_meta("VAR") // default is "?"
             .default_group("main") // default is "global"
             .flag_symbol("F"); // default is "⚑"

  Parser<0> cli("ex01",
                "ex01 description",
                "v0.0.1",
                "J. Doe");

  cli.add_param("-s", "short param");
  cli.add_param("--long", "long param")->def("long_default")->meta("LONG");
  cli.add_param("-b/--both", "short and long")->def("both_default");

  cli.add_param("-f/--flag", "flag param, without value")->as_flag();

  cli.add_common();

  cli.set_positional_bounds(1, 2, "[FILES]...", "Positionals help.");
  
  BCLI_PARSE(cli, argc, argv)

  const std::vector<std::string> pos = cli.get_positionals();

  std::cerr << "-s " << cli.getp("s")->value() << std::endl;
  std::cerr << "--long " << cli.getp("long")->value() << std::endl;
  std::cerr << "-b/--both " << cli.getp("b")->value() << std::endl;
  std::cerr << "-f/--flag " << cli.getp("f")->is_set() << std::endl;

  std::cerr << "Pos: ";
  for (auto& s : pos)
    std::cerr << s << " ";
  std::cerr << std::endl;
}

// ex01 v0.0.1
// 
// DESCRIPTION
//   ex01 description
// 
// USAGE
//   ex01 -s <VAR> [--long <LONG>] [-b/--both <VAR>] [-f/--flag] [-h/--help] [--version] [FILES]...
// 
// POSITIONALS
//   Positionals help.
// 
// OPTIONS
//   [main] - global parameters
//     -s        - short param 
//        --long - long param {long_default}
//     -b --both - short and long {both_default}
//     -f --flag - flag param, without value [F]
// 
//   [common]
//     -h --help    - Show this message and exit. [F]
//        --version - Show version and exit. [F]