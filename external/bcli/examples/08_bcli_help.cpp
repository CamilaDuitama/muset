#include<bcli/bcli.hpp>

using namespace bc;

std::string help()
{
  std::stringstream ss;
  ss << "Custom help" << "\n";
  ss << "-s short param" << "\n";
  ss << "--long long param" << "\n";
  ss << "-b, --both short and long" << "\n";
  ss << "-f, --flag no param" << "\n";
  return ss.str();
}

int main(int argc, char* argv[])
{
  Parser<0> cli("ex08",
                "ex08 description",
                "v0.0.1",
                "J. Doe");

  cli.add_param("-s", "short param");
  cli.add_param("--long", "long param")->def("long_default")->meta("LONG");
  cli.add_param("-b/--both", "short and long")->def("both_default");

  cli.add_param("-f/--flag", "flag param, without value")->as_flag();

  cli.add_common();

  cli.set_help(&help);
  BCLI_PARSE(cli, argc, argv)
}

// Custom help
// -s short param
// --long long param
// -b, --both short and long
// -f, --flag no param