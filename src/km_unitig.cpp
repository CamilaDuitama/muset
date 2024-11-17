#include <algorithm>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <filesystem>

#include <kseq++/seqio.hpp>
#include "../external/sshash/dictionary.hpp"

#include "common.h"


int main_unitig(int argc, char **argv) {

  std::size_t ksize = 31;
  std::size_t msize = 15;
  std::size_t nb_threads = 1;
  std::string out_fname;
  bool out_writeseq = false;
  bool help_opt = false;

  int c;
  while ((c = getopt(argc, argv, "k:m:o:t:sh")) != -1) {
    switch (c) {
      case 'k':
        ksize = std::strtoul(optarg, NULL, 10);
        break;
      case 'm':
        msize = std::strtoul(optarg, NULL, 10);
        break;
      case 'o':
        out_fname = optarg;
        break;
      case 't':
        nb_threads = std::max((long)1, std::strtol(optarg, NULL, 10));
        break;
      case 's':
        out_writeseq = true;
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

  if(argc-optind != 2 || help_opt) {
    std::cout << "Usage: kmat_tools unitig [options] <unitigs.fasta> <kmer_matrix>\n\n";
    std::cout << "Creates a unitig matrix.\n\n";
    std::cout << "Options:\n";
    std::cout << "  -k INT   k-mer size (must be <= 63) [31]\n";
    std::cout << "  -m INT   minimizer length (must be < k) [15]\n";
    std::cout << "  -o FILE  write unitig matrix to FILE [stdout]\n";
    std::cout << "  -t INT   number of threads [1]\n";
    std::cout << "  -s       write the unitig sequence as first column instead of the identifier\n";
    std::cout << "  -h       print this help message\n";
    return 0;
  }

  std::string utg_file = argv[optind];
  if(!std::filesystem::exists(utg_file.c_str())) {
    std::cerr << "[error] unitig file \"" << utg_file << "\" does not exist" << std::endl;
    return 1;
  }

  std::string mat_file = argv[optind+1];
  if(!std::filesystem::exists(mat_file.c_str())) {
    std::cerr << "[error] matrix file \"" << mat_file << "\" does not exist" << std::endl;
    return 1;
  }

  if(ksize <= 0) { 
    std::cerr << "[error] -k parameter must be greater than zero" << std::endl;
    return 1;
  } else if (ksize > 63) {
    std::cerr << "[error] -k parameter must be at most 63" << std::endl;
    return 1;
  }

  if(msize <= 0) { 
    std::cerr << "[error] -m parameter must be greater than zero" << std::endl;
    return 1;
  } else if (msize >= ksize) {
    std::cerr << "[error] -m parameter must be smaller than k-mer size" << std::endl;
    return 1;
  }

  // print some info to stderr

  std::cerr << "[info] k-mer length: " << ksize << std::endl;
  std::cerr << "[info] minimizer length: " << msize << std::endl;
  std::cerr << "[info] threads: " << nb_threads << std::endl;

  // build sshash-based dictionary of k-mers

  std::cerr << "[info] building k-mer dictionary"  << std::endl;

  sshash::dictionary kmer_dict;
  {
    // std::ofstream ofs("sshash.log", std::ios::out);
    // std::streambuf *coutbuf = std::cout.rdbuf();
    // if (ofs.good()) {
    //   std::cout.rdbuf(ofs.rdbuf());
    // }

    sshash::build_configuration build_config;
    build_config.k = ksize;
    build_config.m = msize;
    build_config.c = 5.0;
    build_config.pthash_threads = nb_threads;
    build_config.canonical_parsing = true;
    build_config.verbose = false;

    kmer_dict.build(utg_file, build_config);

    // std::cout.rdbuf(coutbuf);
  }

  std::cerr << "[info] unitigs processed: " << kmer_dict.num_contigs() << std::endl;
  std::cerr << "[info] k-mers processed: " << kmer_dict.size() << std::endl;

  // process matrix file

  FILE *mat = fopen(mat_file.c_str(),"r");
  if(mat == NULL) {
    std::cerr << "[error] cannot open matrix file \"" << mat_file <<"\"\n";
    return 1;
  }

  char *kmer = (char *)calloc(ksize+1,1);
  char *line = NULL;
  size_t line_size = 0;
  size_t line_count = 0;

  bool has_kmer = next_kmer_and_line(kmer, ksize, &line, &line_size, mat);
  std::size_t n_samples = has_kmer ? samples_number(line) : 0;
  fprintf(stderr,"[info] samples: %lu\n", n_samples);

  using sample_t = std::pair<uint32_t,uint32_t>;
  std::vector<std::vector<sample_t>> utg_samples(kmer_dict.num_contigs());
  for(auto& counts: utg_samples) {
    counts.resize(n_samples);
  }

  while(has_kmer) {
    line_count++;

    auto res = kmer_dict.lookup_advanced(kmer);
    if (res.kmer_id == sshash::constants::invalid_uint64) {
        has_kmer = next_kmer_and_line(kmer, ksize, &line, &line_size, mat);
        continue;
    }

    std::size_t utg_id = res.contig_id;
    auto& counts = utg_samples[utg_id];
    char *start = second_column(line);
    int c = 0; char *tok = strtok(start," \t\n");
    while(tok) {
      uint32_t num = strtoul(tok,NULL,10);
      counts[c].first = add_sat(counts[c].first, uint32_t{num > 0});
      counts[c].second = add_sat(counts[c].second, num);
      tok = strtok(NULL," \t\n");
      c++;
    }

    has_kmer = next_kmer_and_line(kmer, ksize, &line, &line_size, mat);
  }

  free(kmer);
  free(line);
  fclose(mat);

  std::ostream* fpout = &std::cout;
  std::ofstream ofs;
  if(!out_fname.empty()) {
    ofs.open(out_fname.c_str());
    if(!ofs.good()) {
        std::cerr << "[error] cannot open output file \"" << out_fname << "\"\n";
        return 1;
    }
    fpout = &ofs;
  }

  *fpout << std::fixed << std::setprecision(2);

  // write output
  std::cerr << "[info] writing unitig matrix"  << std::endl;

  klibpp::KSeq unitig;
  klibpp::SeqStreamIn utg_ssi(utg_file.c_str());

  for(uint64_t utg_id=0; utg_ssi >> unitig; utg_id++) {  

    *fpout << (out_writeseq ? unitig.seq : unitig.name);
    
    auto& counts = utg_samples[utg_id];
    
    std::size_t utg_nb_kmers = unitig.seq.length()-ksize+1;
    for(auto p : counts) {
        double frac = (1.0 * p.first)/utg_nb_kmers;
        double avg_coverage = (1.0 * p.second)/utg_nb_kmers;
        *fpout << ' ' << avg_coverage << ';' << frac;
    }

    *fpout << '\n';
  }

  if(!out_fname.empty()) {
    ofs.close();
  }

  return 0;
}
