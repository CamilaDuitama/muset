#include <bcli/bcli.hpp>

int main(int argc, char* argv[])
{
  bc::Parser<0> cli("kmtricks",
                    "kmtricks cli",
                     "v0.0.3",
                     "Téo Lemane");


  cli.add_param("-f/--file", "fof that contains path of read files")
    ->checker(bc::check::is_file)->meta("FILE");
  
  cli.add_param("-d/--run-dir", "runtime directory")->meta("DIR");

  cli.add_param("-k/--kmer-size", "size of k-mers")
    ->checker(bc::check::is_number)->def("31")->meta("INT");
  
  cli.add_param("--count-abundance-min", "min abundance for solid k-mers")
    ->checker(bc::check::is_number)->def("2")->meta("INT");
  
  cli.add_param("--abundance-max", "max abundance for solid k-mers")
    ->checker(bc::check::is_number)->def("3000000")->meta("INT");
  
  cli.add_param("--max-count", "max count")
    ->checker(bc::check::is_number)->def("255")->meta("INT");
  
  cli.add_param("--max-memory", "max memory per core (in mb)")
    ->checker(bc::check::is_number)->def("8000")->meta("INT");

  bc::param_t mode = cli.add_param("-m/--mode", "output matrix format: [bin|ascii|pa|bf|bf_trp]")
    ->checker(bc::check::f::in("bin|ascii|pa|bf|bf_trp"))->def("bin")->meta("STR");

  cli.add_param("--log-files", "log file: [repart|superk|count|merge|split]")
    ->checker(bc::check::f::in("repart|superk|count|merge|split"))->def("")->meta("STR");

  cli.add_param("--nb-cores", "number of cores")
    ->checker(bc::check::is_number)->def("8")->meta("INT");

  cli.add_param("--keep-tmp", "keep tmp files")->as_flag();
  cli.add_param("--lz4", "compress tmp files")->as_flag();
  
  
  cli.add_group("merge", "merge options, see README");
  
  cli.add_param("--merge-abundance-min", "min abundance for solid k-mers")
    ->checker(bc::check::is_number)->def("1")->meta("INT");
  
  cli.add_param("--recurrence-min", "min recurrence threshold for solid k-mers")
    ->checker(bc::check::is_number)->def("1")->meta("INT");
  
  cli.add_param("--save-if", "keep a non-solid kmer if it's solid in X dataset")
    ->checker(bc::check::is_number)->def("0")->meta("INT");

  cli.add_param("--skip-merge", "skip merge, only with --mode bf")
    ->as_flag()->depends_on(bc::check::f::in(bc::FLAG_VALUE), mode, bc::check::f::in("bf"));


  cli.add_group("pipeline control", {});

  cli.add_param("--until", "run until step: [repart|superk|count|merge|split]")
    ->checker(bc::check::f::in("repart|superk|count|merge|split|all"))
    ->def("all")->meta("STR");
  
  cli.add_param("--only", "run only step: [repart|superk|count|merge|split]")
    ->checker(bc::check::f::in("repart|superk|count|merge|split|all"))
    ->def("all")->meta("STR");

  cli.add_group("advanced performance tweaks", {});
  
  cli.add_param("--minimizer-type", "minimizer type (0=lexi, 1=freq)")
    ->checker(bc::check::f::range(0, 1))->def("0")->meta("INT");

  cli.add_param("--minimizer-size", "size of minimizer")
    ->checker(bc::check::is_number)->def("10")->meta("INT");

  cli.add_param("--repartition-type", "minimizer repartition (0=unordered, 1=ordered)")
    ->checker(bc::check::f::range(0, 1))->def("0")->meta("INT");
  
  cli.add_param("--nb-partitions", "number of partitions (0=auto)")
    ->checker(bc::check::is_number)->def("0")->meta("INT");

  cli.add_group("hash mode configuration", {});
  
  cli.add_param("--hasher", "sabuhash | xor")
    ->checker(bc::check::f::in("sabuhash|xor"))->def("xor")->meta("STR");
  
  cli.add_param("--split", "split matrix, only with --mode bf_trp [sdsl|howde]")
    ->checker(bc::check::f::in("sdsl|howde|none"))->def("none")->meta("STR")
    ->depends_on(bc::check::f::in("sdsl|howde"), mode, bc::check::f::in("bf_trp"));
  
  cli.add_common();

  BCLI_PARSE(cli, argc, argv)
}

/*

kmtricks v0.0.3

DESCRIPTION
  kmtricks cli

USAGE
  kmtricks -f/--file <FILE> -d/--run-dir <DIR> [-k/--kmer-size <INT>] [--count-abundance-min <INT>] 
           [--abundance-max <INT>] [--max-count <INT>] [--max-memory <INT>] 
           [-m/--mode <STR>] [--log-files <STR>] [--nb-cores <INT>] 
           [--merge-abundance-min <INT>] [--recurrence-min <INT>] [--save-if <INT>] 
           [--until <STR>] [--only <STR>] [--minimizer-type <INT>] [--minimizer-size <INT>] 
           [--repartition-type <INT>] [--nb-partitions <INT>] [--hasher <STR>] 
           [--split <STR>] [--keep-tmp] [--lz4] [--skip-merge] [-h/--help] [-v/--verbose] 
           [-d/--debug] [--version] 

OPTIONS
  [global] - global parameters
    -f --file                - fof that contains path of read files 
    -d --run-dir             - runtime directory 
    -k --kmer-size           - size of k-mers {31}
       --count-abundance-min - min abundance for solid k-mers {2}
       --abundance-max       - max abundance for solid k-mers {3000000}
       --max-count           - max count {255}
       --max-memory          - max memory per core (in mb) {8000}
    -m --mode                - output matrix format: [bin|ascii|pa|bf|bf_trp] {bin}
       --log-files           - log file: [repart|superk|count|merge|split] 
       --nb-cores            - number of cores {8}
       --keep-tmp            - keep tmp files [⚑]
       --lz4                 - compress tmp files [⚑]

  [merge] - merge options, see README
     --merge-abundance-min - min abundance for solid k-mers {1}
     --recurrence-min      - min recurrence threshold for solid k-mers {1}
     --save-if             - keep a non-solid kmer if it's solid in X dataset {0}
     --skip-merge          - skip merge, only with --mode bf [⚑]

  [pipeline control]
     --until - run until step: [repart|superk|count|merge|split] {all}
     --only  - run only step: [repart|superk|count|merge|split] {all}

  [advanced performance tweaks]
     --minimizer-type   - minimizer type (0=lexi, 1=freq) {0}
     --minimizer-size   - size of minimizer {10}
     --repartition-type - minimizer repartition (0=unordered, 1=ordered) {0}
     --nb-partitions    - number of partitions (0=auto) {0}

  [hash mode configuration]
     --hasher - sabuhash | xor {xor}
     --split  - split matrix, only with --mode bf_trp [sdsl|howde] {none}

  [common]
    -h --help    - Show this message and exit. [⚑]
    -v --verbose - Verbose mode. [⚑]
    -d --debug   - Debug mode. [⚑]
       --version - Show version and exit. [⚑]

*/