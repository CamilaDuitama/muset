#pragma once
#include <string>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <filesystem>
#include <iomanip>
#include <unistd.h>

namespace fs = std::filesystem;
namespace cfr
{
const std::string cfr_version = CFRCAT_VERSION;
const std::string cfr_str_version = "cfrcat v" + cfr_version;

void show_help()
{
  std::cerr << "Concatenate FILE(s) to OUT_FILE." << std::endl;
  std::cerr << "Options are applied on each file." << std::endl;
  std::cerr << std::endl;
  std::cerr << "usage:" << std::endl;
  std::cerr << "   cfrcat -o OUT_FILE [OPTIONS]... [FILES]..." << std::endl;
  std::cerr << std::endl;
  std::cerr << "options:" << std::endl;
  std::cerr << "   -o OUT_FILE - output path (@ for stdout but needs redirection," << std::endl;
  std::cerr << "                 output in fd referring to a terminal is not supported)." << std::endl;
  std::cerr << "   -r FILE     - a header that will be added to the output." << std::endl;
  std::cerr << "   -l INT      - ignore the l first lines." << std::endl;
  std::cerr << "   -b INT      - ignore the n first bytes." << std::endl;
  std::cerr << "   -c INT      - ignore bytes until c (in [0, 255]) has been seen n times, need -n." << std::endl;
  std::cerr << "   -n INT      - see -c." << std::endl;
  std::cerr << "   -a          - if output file exists, append without erase." << std::endl;
  std::cerr << "   -d          - Display parameters." << std::endl;
  std::cerr << "   -h          - Display this message and exit." << std::endl;
  std::cerr << "   -v          - Display version information and exit." << std::endl;

  std::cerr << std::endl;
  std::cerr << "examples:" << std::endl;
  std::cerr << "   Concat f1 and f2 to output but ignoring two first lines in both:" << std::endl;
  std::cerr << "      cfrcat -o output.txt -l 2 f1.txt f2.txt" << std::endl;
  std::cerr << "   Concat f1 to output but ignoring two first bytes in f1:" << std::endl;
  std::cerr << "      cfrcat -o output.txt -b 2 f1.txt" << std::endl;
  std::cerr << "   Concat f1 and f2 to output but ignoring all bytes before have seen 10 'A' (65) in both:" << std::endl;
  std::cerr << "      cfrcat -o output.txt -c 65 -n 10 f1.txt f2.txt" << std::endl;
}

void die(const std::string& msg, bool help = true)
{
  if (!msg.empty()) std::cerr << msg << std::endl;
  if (help) show_help();
  exit(EXIT_FAILURE);
}

struct cfrcat_params
{
public:
  cfrcat_params()
    : skip_lines(0),
      skip_bytes(0),
      skip_char(-1),
      n(0),
      append(false),
      out(""),
      fof(""),
      header("")
  {}

  void display()
  {
    std::cerr << "skip_lines = " << std::to_string(skip_lines) << std::endl;
    std::cerr << "skip_bytes = " << std::to_string(skip_bytes) << std::endl;
    std::cerr << "skip_char  = " << std::to_string(skip_char) << std::endl;
    std::cerr << "n          = " << std::to_string(n) << std::endl;
    std::cerr << "out        = " << out << std::endl;
    std::cerr << "header     = " << header << std::endl;
    std::cerr << "append     = " << std::boolalpha << append << std::endl;
  }
  
  void is_valid()
  {
    if (out.empty())
      die("ERROR: -o is required.", false);
    
    if (!header.empty())
      if (!fs::exists(header))
        die(header + " doesn't exist.", false);
    
    if (ins.empty())
      die("ERROR: input files are missing.", false);
    else
      for (auto& i: ins)
        if (!fs::exists(i))
          die(i + " doesn't exist.", false);
    
    if ((!n && skip_char != -1) || (n > 0 && skip_char == -1))
      die("-c and -n must be used together.", false);
  }

  int skip_lines;
  int skip_bytes;
  int skip_char;
  int n;
  bool append;
  std::string header;
  std::string out;
  std::string fof;
  std::vector<std::string> ins;
};

typedef struct cfrcat_params cfrcat_params;

cfrcat_params parse(int argc, char* argv[])
{
  if (argc < 2) die("");
  
  cfrcat_params p;
  int c;

  std::stringstream msg;
  bool display = false;
  while ((c = getopt(argc, argv, ":o:l:b:c:n:f:r:dahv")) != -1)
  {
    switch (c)
    {
    case 'o': p.out = optarg; break;
    case 'l': p.skip_lines = std::stoi(optarg); break;
    case 'b': p.skip_bytes = std::stoi(optarg); break;
    case 'n': p.n = std::stoi(optarg); break;
    case 'c': p.skip_char = std::stoi(optarg); break;
    case 'f': p.fof = optarg; break;
    case 'r': p.header = optarg; break;
    case 'a': p.append = true; break;
    case 'd': display = true; break;
    case 'h': die("");
    case 'v': die(cfr_str_version, false);
    case ':':
      msg << "ERROR: -" << (char)optopt << " needs a value";
      die(msg.str());
    case '?':
      msg << "ERROR: Unknow option -" << (char)optopt;
      die(msg.str());
    default:
      exit(EXIT_FAILURE);
    }
  }

  msg.clear();
  for (; optind < argc; optind++)
    p.ins.push_back(argv[optind]);

  p.is_valid();

  if (display)
    p.display();

  return p;
}

};