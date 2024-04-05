#!/bin/bash

# Script Name: kmer_2_unitig_matrix
# Author: Camila Duitama, Francesco Andreace, Riccardo Vicedomini
# Date: 2nd of April 2024
# Description: bash script to create unitig matrix from fof.txt file and set of FASTA/FASTQ files using kmat_tools


SCRIPT_NAME=$(basename $0)

#SETTING DEFAULT PARAMETERS
DEFAULT_KSIZE=31
DEFAULT_WSIZE=14
DEFAULT_THREADS=4
DEFAULT_UTG_MIN_SIZE=100
DEFAULT_MIN_COUNT=1
DEFAULT_KMAT_ABUNDANCE=1
DEFAULT_KMAT_MIN_ZEROS=0
DEFAULT_KMAT_MIN_NONZEROS=0
DEFAULT_UTG_MIN=100
DEFAULT_OUTDIR=output

USAGE=$'\nUsage: '"${SCRIPT_NAME}"' [-k KMER-SIZE] [-t NUM-THREADS] [-u MIN-UTG-SIZE] [-c MIN-COUNT] [-o OUT-DIR] <input_seqfile> <name>

Arguments:
     -h              print this help and exit
     -k              kmer size for ggcat and kmtricks (default:31)
     -t              number of cores (default:4)
     -u              minimum size of the unitigs to be retained in the final matrix (default:100)
     -c              minimum count of kmers to be retained (default:1)
     -a              minimum abundance for kmtool (default:1)
     -n              minimum number of zero columns for each kmers (for kmtools) (default:0)
     -N              minimum number of non-zero columns  for each kmers (for kmtools) (default:0)
     -o              output directory (default:output)


Positional arguments:
     <input_file>         input seqfile (fof)
     <name>               name for output files
'


#PARSING INPUT OPTIONS
k_len=$DEFAULT_KSIZE
thr=$DEFAULT_THREADS
utg_len=$DEFAULT_UTG_MIN_SIZE
min_count=$DEFAULT_MIN_COUNT
output_dir=$DEFAULT_OUTDIR
kmt_abundance=$DEFAULT_KMAT_ABUNDANCE
kmt_min_zeros=$DEFAULT_KMAT_MIN_ZEROS
kmt_min_nonzeros=$DEFAULT_KMAT_MIN_NONZEROS

while getopts ":hktuconNa" opt; do
    case "${opt}" in
	h) $(>&2 echo "${USAGE}")
	   exit 0
	   ;;
	k) k_len=${OPTARG}
	   ;;
   t) thr=${OPTARG}
	   ;;
   u) utg_len=${OPTARG}
	   ;;
	c) maxlen=${OPTARG}
	   ;;
   o) output_dir=${OPTARG}
      ;;
   n) kmt_min_zeros=${OPTARG}
      ;;
   N) kmt_min_nonzeros=${OPTARG}
      ;;
   a) kmt_abundance=${OPTARG}
      ;;

    esac
done

#ADDING INPUT FILE AND OUTPUT FOLDER

input_file=""
#name=""

if [[ $# -lt $((${OPTIND})) ]] # + 1
then
    (>&2 echo "ERROR: Wrong number of arguments.")
    (>&2 echo "")
    (>&2 echo "${USAGE}")
    exit 1
fi

input_file=${@:$OPTIND:1}
#name=${@:$OPTIND+1:1}

start=`date +%s%3N`
#building output directory if not present already
if [ -d $output_dir ]; then rm -rf $output_dir; fi

#Step 1: Produce sorted input k-mer matrix from fof.txt file
runtime=$((`date +%s%3N`-$start)) && echo "[PIPELINE]::[running kmtricks]::[$runtime ms]"

kmtricks pipeline --file $input_file --kmer-size $k_len --hard-min $min_count --mode kmer:count:bin --cpr --run-dir $output_dir

kmtricks aggregate --matrix kmer --format text --cpr-in --sorted --output $output_dir/sorted_matrix.txt --run-dir $output_dir 

# Step 2: Perform some basic filtering with km_basic_filter
runtime=$((`date +%s%3N`-$start)) && echo "[PIPELINE]::[filtering kmers]::[$runtime ms]"

kmtools filter $output_dir/sorted_matrix.txt -a $kmt_abundance -n $kmt_min_zeros -N $kmt_min_nonzeros -o $output_dir/sorted_filtered_matrix.txt

#Step 3: Build unitigs
# Step 3.1: Convert sorted filtered matrix into a fasta file
runtime=$((`date +%s%3N`-$start)) && echo "[PIPELINE]::[exporting filtered kmers]::[$runtime ms]"

kmtools fasta $output_dir/sorted_filtered_matrix.txt -o $output_dir/kmer_matrix.fasta

#Step 3.2: Use ggcat to produce unitigs
runtime=$((`date +%s%3N`-$start)) && echo "[PIPELINE]::[generating unitigs]::[$runtime ms]"

ggcat build $output_dir/kmer_matrix.fasta -o $output_dir/unitigs.fa -j $thr -s $DEFAULT_MIN_COUNT -k $k_len 

#Step 3.3 & 3.4: Filter unitigs that are <100 nt and rename filtered unitigs to be labelled sequentially
awk -v "utg_l=$utg_len" 'BEGIN {count=0} /^>/ { next } { if(length($0) >= utg_l) print ">" count++ "\n" $0 }' $output_dir/unitigs.fa > $output_dir/unitigs_filtered_sequential.fa

#Step 4: Build unitig matrix 
runtime=$((`date +%s%3N`-$start))
echo "[PIPELINE]::[building unitigs matrix]::[$runtime ms]"

kmtools unitig -c $min_count -k $k_len -o $output_dir/unitig_matrix $output_dir/unitigs_filtered_sequential.fa $output_dir/sorted_filtered_matrix.txt

runtime=$((`date +%s%3N`-$start)) && echo "[PIPELINE]::[END]::[$runtime ms]"
exit 0