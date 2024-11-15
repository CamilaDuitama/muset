#include<bcli/bcli.hpp>

using namespace bc;

int main(int argc, char* argv[])
{
  Parser<0> cli("ex03",
                "ex03 description, try: 03_bcli_checker",
                "v0.0.1",
                "J. Doe");

  cli.add_param("-i", "a integer")
    ->checker(check::is_number);

  cli.add_param("-d", "a dir")
    ->def(".")
    ->checker(check::is_dir);

  auto is_localhost = [](const std::string& s, const std::string v) {
    return std::make_tuple(v == "localhost", s + " != localhost");
  };
  
  cli.add_param("-a", "ip")
    ->def("localhost")
    ->checker(is_localhost);


  // Factories (configurable checkers)
  cli.add_param("-s/--size", "a size between 10 and 100")
    ->def("20")
    ->checker(check::f::range(10,100));

  auto seems_jpeg = check::f::ext("jpg|JPEG|jpeg"); // jpeg checker
  auto is_jpeg = check::f::check_magic<4>("jpg", {0xFF, 0xD8, 0xFF, 0xDB}); // jpeg magic checker
  
  // 
  cli.add_param("-j/--jpeg", "a jpeg image")
    ->def("ok.jpeg")
    ->meta("JPEG")
    ->checker(seems_jpeg)
    ->checker(is_jpeg);
  
  cli.add_common();

  BCLI_PARSE(cli, argc, argv)
}

// ex03 v0.0.1
// 
// DESCRIPTION
//   ex03 description, try: 03_bcli_checker
// 
// USAGE
//   ex03 -i <?> [-d <?>] [-a <?>] [-s/--size <?>] [-j/--jpeg <JPEG>] [-h/--help] [-v/--verbose] 
//        [-d/--debug] [--version] 
// 
// OPTIONS
//   [global] - global parameters
//     -i        - a integer 
//     -d        - a dir {.}
//     -a        - ip {localhost}
//     -s --size - a size between 10 and 100 {20}
//     -j --jpeg - a jpeg image {ok.jpeg}
// 
//   [common]
//     -h --help    - Show this message and exit. [⚑]
//     -v --verbose - Verbose mode. [⚑]
//     -d --debug   - Debug mode. [⚑]
//        --version - Show version and exit. [⚑]