#include <bcli/bcli.hpp>

using namespace bc;

int main(int argc, char* argv[])
{
  Parser<0> cli("ex04",
                "ex04 desc",
                "v0.0.1",
                "J. Doe");


  param_t p = cli.add_param("-p/--param", "param help");
  
  // ok: 04_bcli_deps_bans -p 110 -d 8
  // throw: 04_bcli_deps_bans -p 110 -d 20
  // ok: 04_bcli_deps_bans -p 125 -d 20
  cli.add_param("-d/--dep", "dep help")->depends_on(
    check::f::range(10, 100), p, check::f::range(120, 140)
  )->def("120");

  auto c = [](const std::string& s, const std::string& v) {
    return std::make_tuple(v == "15", s + " should not be equal to 15.");
  };

  // ok: 04_bcli_deps_bans -p 15 -b 9
  // throw: 04_bcli_deps_bans -p 15 -b 20
  cli.add_param("-b/--ban", "ban help")->banned(
    check::f::range(10, 100), p, c
  )->def("9");

  cli.add_common();

  BCLI_PARSE(cli, argc, argv)

  return 0;
}

// ex04 v0.0.1
//
// DESCRIPTION
//   ex04 desc
//
// USAGE
//   ex04 -p/--param <?> [-d/--dep <?>] [-b/--ban <?>] [-h/--help] [-v/--verbose] [-d/--debug] 
//        [--version] 
//
// OPTIONS
//   [global] - global parameters
//     -p --param - param help 
//     -d --dep   - dep help {120}
//     -b --ban   - ban help {9}
//
//   [common]
//     -h --help    - Show this message and exit. [⚑]
//     -v --verbose - Verbose mode. [⚑]
//     -d --debug   - Debug mode. [⚑]
//        --version - Show version and exit. [⚑]