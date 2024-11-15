#include<bcli/bcli.hpp>

using namespace bc;

int main(int argc, char* argv[])
{
  Parser<0> cli("ex02",
                "ex02 description, try: ./05_bcli_settter -s 10",
                "v0.0.1",
                "J. Doe");

  std::string s;
  cli.add_param("-s", "set string")->def("10")->setter(s);
  
  int i;
  cli.add_param("-i", "set int")->def("32")->setter(i);

  float f;
  cli.add_param("-f", "set float")->def("18.6")->setter(f);
  
  std::vector<std::string> vec;
  auto my_setter = [&vec] (const std::string& v) {
    vec = utils::split(v, ',');
  };

  cli.add_param("-v", "set vector")->def("1,2,3")->setter_c(my_setter);

  cli.add_common();

  BCLI_PARSE(cli, argc, argv)

  std::cerr << s << std::endl;
  std::cerr << i << std::endl;
  std::cerr << f << std::endl;

  for (auto& e: vec)
    std::cerr << e << " ";
  std::cerr << std::endl;
}

// ex02 v0.0.1
//
// DESCRIPTION
//   ex02 description, try: ./05_bcli_settter -s 10
//
// USAGE
//   ex02 [-s <?>] [-i <?>] [-f <?>] [-v <?>] [-h/--help] [-v/--verbose] [-d/--debug] [--version] 
//
// OPTIONS
//   [global] - global parameters
//     -s  - set string {10}
//     -i  - set int {32}
//     -f  - set float {18.6}
//     -v  - set vector {1,2,3}
//
//   [common]
//     -h --help    - Show this message and exit. [⚑]
//     -v --verbose - Verbose mode. [⚑]
//     -d --debug   - Debug mode. [⚑]
//        --version - Show version and exit. [⚑]
