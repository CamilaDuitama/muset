#include "common.h"


int main_merge(int argc, char **argv) {

  int ksize = 31;
  char *out_fname = NULL;
  bool use_ktcmp = false, help_opt = false;

  int c;
  while ((c = getopt(argc, argv, "k:o:zh")) != -1) {
    switch (c) {
      case 'k':
        ksize = strtol(optarg, NULL, 10);
        break;
      case 'o':
        out_fname = optarg;
        break;
      case 'z':
        use_ktcmp = true;
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

  if(ksize <= 0) { 
    fprintf(stderr, "Invalid value of k: %d\n",ksize);
    return 1;
  }

  if(argc-optind != 2 || help_opt) {
    fprintf(stdout, "Usage: kmat_tools merge [options] <matrix_1> <matrix_2>\n\n");
    fprintf(stdout, "Merge two input kmer-sorted matrices.\n\n");
    fprintf(stdout, "Options:\n");
    fprintf(stdout, "  -k INT   size of k-mers of input matrices [31]\n");
    fprintf(stdout, "  -o FILE  write output matrix to FILE [stdout]\n");
    fprintf(stdout, "  -z       use kmtricks order of nucleotides: A<C<T<G\n");
    fprintf(stdout, "  -h       print this help message\n");
    return 0;
  }

  FILE *mat_1 = fopen(argv[optind],"r");
  if(mat_1 == NULL) {
    fprintf(stderr,"Cannot open file \"%s\"\n",argv[optind]);
    return 1;
  }

  FILE *mat_2 = fopen(argv[optind+1],"r");
  if(mat_2 == NULL) {
    fprintf(stderr,"Cannot open file \"%s\"\n",argv[optind+1]);
    fclose(mat_1); 
    return 1;
  }

  FILE *outfile = out_fname ? fopen(out_fname,"w") : stdout;
  if(outfile != stdout && outfile == NULL) {
    fprintf(stderr,"Cannot open output file \"%s\"\n",out_fname);
    fclose(mat_1);
    fclose(mat_2);
    return 1;
  }

  char *kmer_1 = (char *)calloc(ksize+1,1);
  char *kmer_2 = (char *)calloc(ksize+1,1);
  char *line_1 = NULL, *line_2 = NULL;
  size_t line_1_size = 0, line_2_size = 0;

  bool has_kmer_1 = next_kmer_and_line(kmer_1, ksize, &line_1, &line_1_size, mat_1);
  size_t n_sample_1 = has_kmer_1 ? samples_number(line_1) : 0;
  fprintf(stderr,"[info] samples in 1st matrix: %lu\n", n_sample_1);

  bool has_kmer_2 = next_kmer_and_line(kmer_2, ksize, &line_2, &line_2_size, mat_2);
  size_t n_sample_2 = has_kmer_2 ? samples_number(line_2) : 0;
  fprintf(stderr,"[info] samples in 2nd matrix: %lu\n", n_sample_2);

  while(has_kmer_1 && has_kmer_2){
    int ret_cmp = use_ktcmp ? ktcmp(kmer_1,kmer_2) : strcmp(kmer_1,kmer_2);
    if(ret_cmp == 0) {
      fputs(kmer_1,outfile);
      fputc(' ',outfile);
      fputs(second_column(line_1),outfile);
      fputc(' ',outfile);
      fputs(second_column(line_2),outfile);
      has_kmer_1 = next_kmer_and_line(kmer_1, ksize, &line_1, &line_1_size, mat_1);
      has_kmer_2 = next_kmer_and_line(kmer_2, ksize, &line_2, &line_2_size, mat_2);
    } else if(ret_cmp < 0) {
      fputs(kmer_1,outfile);
      fputc(' ',outfile);
      fputs(second_column(line_1),outfile);
      for(size_t i=0; i<n_sample_2; ++i){ fputs(" 0",outfile); }
      has_kmer_1 = next_kmer_and_line(kmer_1, ksize, &line_1, &line_1_size, mat_1);
    } else { // ret_cmp > 0
      fputs(kmer_2,outfile);
      for(size_t i=0; i<n_sample_1; ++i){ fputs(" 0",outfile); }
      fputc(' ',outfile);
      fputs(second_column(line_2),outfile);
      has_kmer_2 = next_kmer_and_line(kmer_2, ksize, &line_2, &line_2_size, mat_2);
    }
    fputc('\n',outfile);
  }

  while(has_kmer_1) {
    fputs(kmer_1,outfile);
    fputc(' ',outfile);
    fputs(second_column(line_1),outfile);
    for(size_t i=0; i<n_sample_2; ++i){ fputs(" 0",outfile); }
    fputc('\n',outfile);
    has_kmer_1 = next_kmer_and_line(kmer_1, ksize, &line_1, &line_1_size, mat_1);
  }

  while(has_kmer_2) {
    fputs(kmer_2,outfile);
    for(size_t i=0; i<n_sample_1; ++i){ fputs(" 0",outfile); }
    fputc(' ',outfile);
    fputs(second_column(line_2),outfile);
    fputc('\n',outfile);
    has_kmer_2 = next_kmer_and_line(kmer_2, ksize, &line_2, &line_2_size, mat_2);
  }

  free(kmer_1);
  free(kmer_2);
  free(line_1);
  free(line_2);
  fclose(mat_1);
  fclose(mat_2);
  if(outfile != stdout){ fclose(outfile); }

  return 0;
}
