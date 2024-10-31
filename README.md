<div style="display: flex; align-items: center;">
  <img src="logo.png" alt="Logo" width="200" style="margin: 0;">
</div>

## A pipeline for building an abundance unitig matrix.
MUSET is a software for generating an abundance unitig matrix from a collection of input samples (in FASTA/Q format).
It additionally provides a comprehensive suite of tools (called `kmat tools`) for manipulating k-mer matrices and a script  
for generating a presence-absence unitig matrix.

+ [Installation](#installation)
  - [Conda installation](#conda-installation)
  - [Build from source](#build-from-source)
  - [Build a Singularity image](#build-a-singularity-image)
+ [Usage](#usage)
  - [Input data](#input-data)
    - [I do not have a k-mer matrix](#i-do-not-have-a-k-mer-matrix)
    - [I already have a k-mer matrix](#i-already-have-a-k-mer-matrix)
  - [Output data](#output-data)
  - [k-mer matrix operations](#k-mer-matrix-operations)
  - [I just want a presence-absence unitig matrix](#i-just-want-a-presence-absence-unitig-matrix)
    - [Input data](#input-file)
    - [Output file](#output-file)
+ [Acknowledgements](#acknowledgements)

## Installation

### Conda installation

The simplest way to install `muset` is by creating a conda environment (e.g., `muset_env`) as follows:
```
conda create -n muset_env -c conda-forge bioconda::ggcat camiladuitama::muset
```

To run `muset` remember to activate the corresponding conda environment with:
```
conda activate muset_env
```

You can check if `muset` is correctly installed as follows:
```
cd test
muset fof.txt
```

### Build from source

Requirements:
  - a recent version of GCC (or clang) that supports the C++17 standard
  - cmake >= 3.15
  - [GGCAT](https://github.com/algbio/ggcat?tab=readme-ov-file#installation)

To clone the repository:
```
git clone --recursive https://github.com/camiladuitama/muset.git
```

To build the tool:
```
cd muset
mkdir build && cd build
cmake ..
make
```
Executables will be made available in the `bin` sub-directory relative to the root folder of the repository.

To make the `muset` command available, remember to include the absolute path of MUSET's executables in your PATH environment variable, e.g., adding the following line to your `~/.bashrc` file:
```
export PATH=/absolute/path/to/muset/bin:${PATH}
```

You can check if `muset` is correctly installed as follows:
```
cd test
muset fof.txt
```

### Build a Singularity image

Requirements:
  - Singularity installed on your system. Refer to the [Singularity Installation Guide](https://sylabs.io/guides/3.0/user-guide/installation.html) for detailed instructions.


To build a singularity image (e.g., `muset.sif`):
```
git clone --recursive https://github.com/CamilaDuitama/muset.git
cd muset/singularity
sudo singularity build muset.sif Singularity.def
```

To run `muset` and see the help message, use the following command:
```
singularity exec /path/to/muset.sif muset -h
```

To try `muset` with example data, `cd` to the `test` directory within the repository, then run:
```
singularity exec /path/to/muset.sif muset fof.txt
```

## Usage

````
muset v0.4

DESCRIPTION:
   muset - a pipeline for building an abundance unitig matrix from a list of FASTA/FASTQ files.

USAGE:
   muset [options] INPUT_FILE

OPTIONS:
   -i PATH    skip matrix construction and run the pipeline with a previosuly computed matrix
   -k INT     k-mer size (default: 31)
   -a INT     min abundance to keep a k-mer (default: 2)
   -l INT     minimum size of the unitigs to be retained in the final matrix (default: 2k-1)
   -o PATH    output directory (default: output)
   -m INT     minimizer length  (default: 15)
   -n INT     minimum number of samples from which a k-mer should be absent (mutually exclusive with -f)
   -f FLOAT   fraction of samples from which a k-mer should be absent (default: 0.1, mutually exclusive with -n)
   -N INT     minimum number of samples in which a k-mer should be present (mutually exclusive with -F)
   -F FLOAT   fraction of samples in which a k-mer should be present (default: 0.1, mutually exclusive with -N)
   -t INT     number of cores (default: 4)
   -s         write the unitig sequence in the first column of the output matrix instead of the identifier
   -h         show this help message and exit
   -V         show version number and exit

POSITIONAL ARGUMENTS:
    INPUT_FILE   Input file (fof) containing the description of input samples.
                 It is ignored if -i option is used.

NOTES:
   Options -n and -f are mutually exclusive, as well as options -N and -F.
   When either -n or -f is used, -N or -F must also be provided, and vice versa.
   If none of the -n, -N, -f, -F options are used the last two options are used with their default values.
````

### Input data

### I do not have a k-mer matrix

If you do not have a k-mer matrix ready, make sure to create a "fof" file, that is a file which contains one line per sample with the following syntax:
  - `<Sample ID> : <1.fastq.gz> ; ... ; <N.fastq.gz>`

Files could be in either FASTA or FASTQ format, gzipped or not.
Multiple files per sample can be provided by separating them with a semicolon.

<ins>Example:</ins>
```
A1 : /path/to/fastq_A1_1
B1 : /path/to/fastq_B1_1 ; /with/mutiple/fasta_B1_2
```

You can generate such an input file from a folder containing many input files as follows:
```
ls -1 folder/*  | sort -n -t/ -k2 | xargs -L1 readlink -f | awk '{ print ++count" : "$1 }' >fof.txt
```

Then simply run:
```
muset fof.txt
```

### I already have a k-mer matrix
If you are familiar with `kmtricks` and/or have already produced a k-mer matrix on your own, you can run `muset` with the `-i` option and provide your own input matrix (and skip the possibly long matrix construction).

Make sure to provide a matrix in text format. You can easily output one from a kmtricks run using the command `kmtricks aggregate` with parameters `--matrix kmer --format text`.
By default, `kmtricks` will write it on stdout, so you might want to set the `--output` parameter.
Ex: `kmtricks aggregate --matrix kmer --format text --cpr-in --sorted --output sorted_matrix.txt --run-dir kmtricks_output_dir`

The pipeline can be then run as follows:
```
muset -i sorted_matrix.txt <input_fof.txt>
```
### Output data

The output data of `muset` is a folder with intermediate results and a `unitigs.mat` file, which is an abundance unitig matrix. Each row corresponds to a unitig and each column to a sample. Each entry of the matrix indicates the average abundance and fraction of the unitig k-mers belonging to the sample (separated by a semicolon) Ex:

| Unitig ID | Sample 1 | Sample 2 | Sample 3      | Sample 4      | Sample 5      |
|-----|----------|----------|---------------|---------------|---------------|
| 0   | 0.00;0.00| 0.00;0.00| 0.00;0.00     | 0.00;0.00     | 2.00;1.00     |
| 1   | 2.00;1.00| 2.00;1.00| 2.00;1.00     | 2.00;1.00     | 0.00;0.00     |
| 2   | 0.00;0.00| 0.00;0.00| 0.00;0.00     | 0.00;0.00     | 2.00;1.00     |
| 3   | 0.00;0.00| 0.00;0.00| 0.00;0.00     | 0.00;0.00     | 2.00;1.00     |
| 4   | 2.00;1.00| 2.00;1.00| 2.00;1.00     | 2.00;1.00     | 0.00;0.00     |

**Note:** If instead of the unitig identifier you prefer to have the unitig sequence, run `muset` with the flag `-s`

The average abundance of a unitig $u$ with respect to a sample $S$ (number on the left of the semicolon) is defined as:

$$ A(u, S) = \frac{\sum\limits_{i=1}^{N}{c_i}}{N} $$

where $N$ is the number of k-mers in $u$, and $c_i$ is the abundance of the $i$-th k-mer of $u$ in sample $S$.

The fraction of k-mers in a unitig $u$ that are present in a sample $S$ (number on the right of the semicolon) is defined as:

$$ f(u, S) = \frac{\sum\limits_{i=1}^{N}{x_i}}{N} $$

where $N$ is the number of k-mers in $u$, and $x_i$ is a binary variable that is 1 when the $i$-th k-mer is present in sample $S$ and 0 otherwise.


### K-mer matrix operations

MUSET includes a `kmat_tools`, an auxiliary executable allowing to perform some basic operations on a (text) k-mer matrix.

```
kmat_tools v0.4

DESCRIPTION
  kmat_tools - a collection of tools to process text-based k-mer matrices

USAGE
  kmat_tools <command> <arguments>

COMMANDS
  diff     - difference between two sorted k-mer matrices
  fasta    - output a k-mer matrix in FASTA format
  fafmt    - filter a FASTA file by length and write sequences in single lines
  filter   - filter a k-mer matrix by selecting k-mers that are potentially differential
  merge    - merge two input sorted k-mer matrices
  reverse  - reverse complement k-mers in a matrix
  select   - select only a subset of k-mers
  unitig   - build a unitig matrix
  version  - print version
```

### I just want a presence-absence unitig matrix
MUSET includes also `muset_pa`, an auxiliary executable that generates a presence-absence unitig matrix in text format from a list of input samples using ggcat and kmat_tools.

```
muset_pa v0.4

DESCRIPTION:
   muset_pa - a pipeline for building a presence-absence unitig matrix from a list of FASTA/FASTQ files.
              this pipeline has fewer parameters than muset and less filtering options as it does not build
              nor use an intermediate k-mer abundance matrix.
              If you wish a 0/1 binary matrix instead of the fraction of kmers from the sample present in the
              unitig, please use the option -r and a value x, 0 < x <=1 as minimum treshold to count a sample
              as present (1).

USAGE:
   muset_pa [options] INPUT_FILE

OPTIONS:
   -k INT     k-mer size (default: 31)
   -a INT     min abundance to keep a k-mer (default: 2)
   -l INT     minimum size of the unitigs to be retained in the final matrix (default: 2k-1)
   -r FLOAT   minimum kmer presence ratio in the unitig for 1/0 
   -o PATH    output directory (default: output)
   -m INT     minimizer length  (default: 15)
   -t INT     number of cores (default: 4)
   -s         write the unitig sequence in the first column of the output matrix instead of the identifier
   -h         show this help message and exit
   -V         show version number and exit

POSITIONAL ARGUMENTS:
    INPUT_FILE   Input file (fof) containing paths of input samples (one per line).
```

#### Input file
The input is a file containing a list of paths (one per line), as required by the `-l` parameter of GGCAT.
Make sure to either specify absolute paths or paths relative to the directory from which you intend to run `muset_pa`.

A simple test example can be run from the `test` directory:
```
cd test
muset_pa -o output_pa fof_pa.txt
```


#### Output file
The pipeline will produce multiple intermediate output files, among which the jsonl dictionary of the colors for each unitig that is normally produced by ggcat.
The pipeline automatically converts it into a unitig matrix in csv format (separated by column). If you choose option -r you will have it in binary format (0/1) else it will report the
percentage of k-mers from each samples inside the unitigs. Samples will have the name of the input files you used.
Here an example  
| Unitig ID | Sample 1 | Sample 2 | Sample 3      | Sample 4      | Sample 5      |
|-----|-----|-----|-----|-----|-----|
| 0   | 0   | 1   |0.23 | 0.3 |  1  |
| 1   | 1   | 1   |  0  | 0.8 | 0.4 |
| 2   | 0.47| 0.2 |  1  |  1  |  0  |
| 3   | 0.8 | 1   |0.78 |  1  | 0.81|
| 4   | 0.79| 1   |  1  | 0.87|0.89 |

In case you use -r 0.8, you will have  
| Unitig ID | Sample 1 | Sample 2 | Sample 3      | Sample 4      | Sample 5      |
|-----|-----|-----|-----|-----|-----|
| 0   |  0  |  1  |  0  |  0  |  1  |
| 1   |  1  |  1  |  0  |  1  |  0  |
| 2   |  0  |  0  |  1  |  1  |  0  |
| 3   |  1  |  1  |  0  |  1  |  1  |
| 4   |  0  |  1  |  1  |  1  |  1  |


## Acknowledgements

MUSET is based on the following libraries (included in the `external` directory along with their license):

- [kseq++](https://github.com/cartoonist/kseqpp): parsing of FASTA file
- [PTHash](https://github.com/jermp/pthash): compact minimal perfect hash
- [SSHash](https://github.com/jermp/sshash): Sparse and Skew Hashing of K-Mers

For building a k-mer matrix and unitigs the following two software are used:

- [kmtricks](https://github.com/tlemane/kmtricks)
- [GGCAT](https://github.com/algbio/ggcat)
