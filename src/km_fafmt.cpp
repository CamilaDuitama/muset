#include <iostream>
#include <fstream>
#include <filesystem>

#include <kseq++/seqio.hpp>

#include "common.h"


int main_fafmt(int argc, char **argv) {

  std::size_t min_length = 0;
  std::string out_file;
  bool help_opt = false;

  int c;
  while ((c = getopt(argc, argv, "l:o:h")) != -1) {
    switch (c) {
      case 'l':
        min_length = strtoul(optarg, NULL, 10);
        break;
      case 'o':
        out_file = optarg;
        break;
      case 'h':
        help_opt = true;
        break;
      case '?':
        return 1;
      default:
        abort();
    }
  }

  if(argc-optind != 1 || help_opt) {
    fprintf(stdout, "Usage: kmat_tools fafmt [options] <input.fa>\n\n");

    fprintf(stdout, "Filter a fasta file by length and writes its sequences in single lines\n\n");

    fprintf(stdout, "Options:\n");
    fprintf(stdout, "  -l INT   minimum sequence length [0]\n");
    fprintf(stdout, "  -o FILE  output FASTA file of k-mers to FILE [stdout]\n");
    fprintf(stdout, "  -h       print this help message\n");
    return 0;
  }

  std::string input_file = argv[optind];
  if(!std::filesystem::exists(input_file.c_str())) {
    std::cerr << "[error] input file \"" << input_file << "\" does not exist" << std::endl;
    return 1;
  }

  std::ostream* fpout = &std::cout;
  std::ofstream ofs;
  if(!out_file.empty()) {
    ofs.open(out_file.c_str());
    if(!ofs.good()) {
        std::cerr << "[error] cannot open output file \"" << out_file << "\"\n";
        return 1;
    }
    fpout = &ofs;
  }

  klibpp::KSeq record;
  klibpp::SeqStreamIn ssi(input_file.c_str());

  while (ssi >> record) {
    if (record.seq.length() < min_length) {
      continue;
    }
    *fpout << '>' << record.name;
    if( !record.comment.empty() ) { *fpout << ' ' << record.comment; }
    *fpout << '\n' << record.seq << '\n';
  }

  if(!out_file.empty()) {
    ofs.close();
  }

  return 0;
}
