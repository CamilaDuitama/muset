#include<bcli/bcli.hpp>

using namespace bc;

std::string mainhelp()
{
  std::stringstream ss;
  ss << "Main help" << "\n";
  ss << "Usage:" << "\n";
  ss << "09_bcli_help_multiple cmd1" << "\n";
  ss << "09_bcli_help_multiple cmd2" << "\n";
  return ss.str();
}

std::string helpcmd1()
{
  std::stringstream ss;
  ss << "Help cmd 1" << "\n";
  return ss.str();
}

std::string helpcmd2()
{
  std::stringstream ss;
  ss << "Help cmd 2" << "\n";
  return ss.str();
}

int main(int argc, char* argv[])
{
  Parser<1> cli("ex09",
                "ex09 description",
                "v0.0.1",
                "J. Doe");

  cli.set_help(&mainhelp);

  cmd_t cmd1 = cli.add_command("cmd1", "cmd1 desc", &helpcmd1);
  cmd1->add_common();
  
  cmd_t cmd2 = cli.add_command("cmd2", "cmd2 desc", &helpcmd2);
  cmd2->add_common();

  BCLI_PARSE(cli, argc, argv)
}