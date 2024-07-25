# Reproducibility of experiments of muset's paper

This repository contains a script to download Sequence Read Archive (SRA) files for `muset` experiments.

## Table of Contents
* [Prerequisites](#prerequisites)
* [Download SRA Accessions](#download-sra-accessions)
* [Running on Small Dataset](#running-on-small-dataset)
* [Running on Large Dataset](#running-on-large-dataset)

## Prerequisites
Before running the script, ensure you have the following tools installed:

- [`gzip`](https://www.gnu.org/software/gzip/gzip.html): A GNU software application used for file compression and decompression.
- `prefetch` and `fasterq-dump` included in the SRA Toolkit, used for downloading and converting SRA data into the FASTQ format. You can find instructions on how to install the SRA Toolkit [here](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit).
- [`muset`](https://github.com/CamilaDuitama/muset/tree/main?tab=readme-ov-file#installation)

## Download SRA Accessions

To download SRA accessions, use the following command:

```
./download_data.sh -i accessions.txt -o output -t 5
```

In this command:

- `-i accessions.txt`: Replace `accessions.txt` with the name of your input file containing the accession list.
- `-o output`: Replace `output` with your preferred output directory.
- `-t 5`: Replace `5` with your preferred number of threads for `fasterq-dump`.

The script generates an output file named `fof.txt`, containing the list of downloaded sequence data files.

## Running on Small Dataset

Once you have the `fof.txt` file that results of running the script `download_accessions.sh`, you can use it to reproduce the results by running the following `kmtricks` commands:

```
kmtricks pipeline --file fof.txt --kmer-size 31 --hard-min 2 --mode kmer:count:bin --cpr --run-dir small_exp -t 20 --recurrence-min 3 --restrict-to-list 1
kmtricks aggregate --matrix kmer --format text --cpr-in --sorted --output small_exp/sorted_matrix.txt --run-dir small_exp -t 20

```

After aggregating the results with kmtricks, you can use the following command to run `muset`:

```
muset -t 20 -l 50 -i small_exp/sorted_matrix.txt
```

## Running on Large Dataset

For large datasets, use the following command to run `muset`':

```
muset -t 20 fof.txt
```


