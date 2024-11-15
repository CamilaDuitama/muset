#include <bcli/bcli.hpp>

int main(int argc, char* argv[])
{
  bc::Parser<1> cli("ex06",
                    "ex06 desc",
                    "v0.0.1",
                    "J. Doe");

  bc::cmd_t merge = cli.add_command("merge", "A merge command");
  
  merge->add_param("-f/--file", "path to file");
  merge->add_common();

  bc::cmd_t split = cli.add_command("split", "A split command");
  split->add_param("-f/--file", "path to file");
  split->add_common();
  
   BCLI_PARSE(cli, argc, argv)
} 

// ex06
// ex06 v0.0.1
// 
// DESCRIPTION
//   ex06 desc
// 
// USAGE
//   ex06 [merge|split]
// 
// COMMANDS
//   merge - A merge command
//   split - A split command

// ex06 merge v0.0.1
// 
// DESCRIPTION
//   A merge command
// 
// USAGE
//   ex06 merge -f/--file <?> [-h/--help] [-v/--verbose] [-d/--debug] [--version] 
// 
// OPTIONS
//   [global]
//     -f --file - path to file 
// 
//   [common]
//     -h --help    - Show this message and exit. [⚑]
//     -v --verbose - Verbose mode. [⚑]
//     -d --debug   - Debug mode. [⚑]
//        --version - Show version and exit. [⚑]
// 

// ex06 split v0.0.1
// 
// DESCRIPTION
//   A split command
// 
// USAGE
//   ex06 split -f/--file <?> [-h/--help] [-v/--verbose] [-d/--debug] [--version] 
// 
// OPTIONS
//   [global]
//     -f --file - path to file 
// 
//   [common]
//     -h --help    - Show this message and exit. [⚑]
//     -v --verbose - Verbose mode. [⚑]
//     -d --debug   - Debug mode. [⚑]
//        --version - Show version and exit. [⚑]