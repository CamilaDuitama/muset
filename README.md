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

The simplest way to install `muset` (and its dependencies) is by creating a conda environment (e.g., `muset_env`) as follows:
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
muset --file fof.txt
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
make -j4
```
Executables will be made available in the `bin` sub-directory of the main repository.

To make the `muset` command available, you might want to include the absolute path of the `bin` directory in your `PATH` environment variable, e.g., adding the following line to your `~/.bashrc` file:
```
export PATH=/absolute/path/to/muset/bin:${PATH}
```

You can check if `muset` is correctly installed as follows:
```
cd test
muset --file fof.txt
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
singularity exec /path/to/muset.sif muset --help
```

To try `muset` with example data, `cd` to the `test` directory within the repository, then run:
```
singularity exec /path/to/muset.sif muset --file fof.txt
```

## Usage

````
muset v0.5

DESCRIPTION
  a pipeline for building an abundance unitig matrix from a list of FASTA/FASTQ files.

USAGE
  muset [--file <FILE>] [-i/--in-matrix <FILE>] [-o/--out-dir <DIR>] [-k/--kmer-size <INT>] 
        [-m/--mini-size <INT>] [-a/--min-abundance <INT>] [-l/--min-unitig-length <INT>] 
        [-r/--min-utg-frac <FLOAT>] [-f/--min-frac-absent <FLOAT>] 
        [-F/--min-frac-present <FLOAT>] [-n/--min-nb-absent <FLOAT>] 
        [-N/--min-nb-present <FLOAT>] [-t/--threads <INT>] [-s/--write-seq] [--out-frac] 
        [-u/--logan] [--keep-temp] [-h/--help] [-v/--version] 

OPTIONS
  [main options]
       --file              - kmtricks-like input file, see README.md. 
    -i --in-matrix         - input matrix (text file or kmtricks directory). 
    -o --out-dir           - output directory. {output}
    -k --kmer-size         - k-mer size. [8, 63]. {31}
    -m --mini-size         - minimizer size. [4, 15]. {15}
    -a --min-abundance     - minimum abundance to keep a k-mer. {2}
    -l --min-unitig-length - minimum unitig length. {2k-1}
    -r --min-utg-frac      - minimum k-mer fraction to set unitig average abundance [0,1]. {0.0}
    -s --write-seq         - write the unitig sequence instead of the identifier in the output matrix [⚑]
       --out-frac          - output an additional matrix containing k-mer fractions. [⚑]
    -u --logan             - input samples consist of Logan unitigs (i.e., with abundance). [⚑]

  [filtering options]
    -f --min-frac-absent  - fraction of samples from which a k-mer should be absent. [0.0, 1.0] {0.1}
    -F --min-frac-present - fraction of samples in which a k-mer should be present. [0.0, 1.0] {0.1}
    -n --min-nb-absent    - minimum number of samples from which a k-mer should be absent (overrides -f). {0}
    -N --min-nb-present   - minimum number of samples in which a k-mer should be present (overrides -F). {0}

  [other options]
       --keep-temp - keep temporary files. [⚑]
    -t --threads   - number of threads. {4}
    -h --help      - show this message and exit. [⚑]
    -v --version   - show version and exit. [⚑]
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
muset --file fof.txt
```

### I already have a k-mer matrix
If you are familiar with [kmtricks](https://github.com/tlemane/kmtricks) and/or have already produced a k-mer matrix on your own, you can run `muset` with the `-i` option and provide your own input matrix (and skip the possibly long matrix construction).
This k-mer matrix can be either a text file (with values separated by a single space) or a [kmtricks](https://github.com/tlemane/kmtricks) output directory.

The pipeline can be then run as follows:
```
muset -i /path/to/input/matrix
```
### Output data

The output of `muset` is a folder with intermediate results and a `unitigs.abundance.mat` file, which is an abundance unitig matrix. Each row corresponds to a unitig and each column to a sample. Each entry of the matrix indicates the average abundance of the unitig k-mers within the corresponding sample Ex:

| Unitig ID | Sample 1 | Sample 2 | Sample 3      | Sample 4      | Sample 5      |
|-----------|----------|----------|---------------|---------------|---------------|
| 0         | 0.00     | 0.00     | 0.00          | 0.00          | 2.00          |
| 1         | 2.00     | 2.00     | 2.00          | 2.00          | 0.00          |
| 2         | 0.00     | 0.00     | 0.00          | 0.00          | 2.00          |
| 3         | 0.00     | 0.00     | 0.00          | 0.00          | 2.00          |
| 4         | 2.00     | 2.00     | 2.00          | 2.00          | 0.00          |

**Note:** 
The sequence of a unitig can be retrieved from the FASTA file `unitigs.fa` stored in the output folder.
If, instead, you prefer to directly have the unitig sequence in the first column, you can run `muset` using the `-s` flag.

The average abundance of a unitig $u$ with respect to a sample $S$ (number on the left of the semicolon) is defined as:

$$ A(u, S) = \frac{\sum\limits_{i=1}^{N}{c_i}}{N} $$

where $N$ is the number of k-mers in $u$, and $c_i$ is the abundance of the $i$-th k-mer of $u$ in sample $S$.

A companion matrix `unitigs.frac.mat`, containing the fraction of the unitig's k-mers belonging to a sample, can be outputted using the `--out-frac` parameter. Moreover the `--min-utg-frac` option allows to output the abundance value in the `unitigs.abundance.mat` file only for a unitig whose fraction is greater than a given threshold (otherwise abundance is set to `0.0`).

The fraction of k-mers in a unitig $u$ that are present in a sample $S$ (number on the right of the semicolon) is defined as:

$$ f(u, S) = \frac{\sum\limits_{i=1}^{N}{x_i}}{N} $$

where $N$ is the number of k-mers in $u$, and $x_i$ is a binary variable that is 1 when the $i$-th k-mer is present in sample $S$ and 0 otherwise.


### K-mer matrix operations

MUSET includes a `kmat_tools`, an auxiliary executable allowing to perform some basic operations on a (text) k-mer matrix.

```
kmat_tools v0.5

DESCRIPTION
  a collection of tools to process text-based k-mer matrices

USAGE
  kmat_tools [convert|diff|fafmt|fasta|filter|merge|reverse|unitig]

COMMANDS
  convert - Convert ggcat jsonl color output into a unitig matrix
  diff    - Difference between two sorted k-mer matrices
  fafmt   - Filter a FASTA file by length and write sequences in single lines
  fasta   - Output a k-mer matrix in FASTA format
  filter  - Filter a matrix by selecting k-mers that are potentially differential.
  merge   - Merge two input text-based kmer-sorted matrices.
  reverse - Reverse-complement k-mers in a k-mer matrix file.
  unitig  - Create a unitig matrix.
```

### I just want a presence-absence unitig matrix
MUSET also includes `muset_pa`, an executable for building a presence-absence unitig matrix in text format from a list of input samples using `ggcat` and `kmat_tools`.

```
muset_pa v0.5

DESCRIPTION
  a pipeline for building a presence-absence unitig matrix from a list of FASTA/FASTQ files.

USAGE
  muset_pa [--file <FILE>] [-o/--out-dir <DIR>] [-k/--kmer-size <INT>] [-m/--mini-size <INT>] 
           [-a/--min-abundance <INT>] [-l/--min-unitig-length <INT>] 
           [-r/--min-utg-frac <FLOAT>] [-t/--threads <INT>] [-s/--write-seq] [-h/--help] 
           [-v/--version] 

OPTIONS
  [main options]
       --file              - list of FASTA/Q files, one per line (see README.md). 
    -o --out-dir           - output directory. {output_pa}
    -k --kmer-size         - k-mer size. [8, 63]. {31}
    -m --mini-size         - minimizer size. [4, 15]. {15}
    -a --min-abundance     - minimum abundance required to keep a kmer. {2}
    -l --min-unitig-length - minimum length required to keep a unitig. {2k-1}
    -r --min-utg-frac      - output a binary matrix; sets a unitig as present (1) when fraction is greater than this threshold [0,1]. {0.8}
    -s --write-seq         - write the unitig sequence instead of the identifier in the output matrix [⚑]

  [other options]
    -t --threads - number of threads. {4}
    -h --help    - show this message and exit. [⚑]
    -v --version - show version and exit. [⚑]
```

#### Input file
The input is a file containing a list of paths (one per line), as required by the `-l` parameter of [GGCAT](https://github.com/algbio/ggcat). Make sure to either specify absolute paths or paths relative to the directory from which you intend to run `muset_pa`.

A simple test example can be run from the `test` directory:
```
cd test
muset_pa -o output_pa --file fof_pa.txt
```


#### Output file

The pipeline will produce several intermediate output files, among which the jsonl dictionary of the colors for each unitig that is normally produced by ggcat. The pipeline automatically converts it into a unitig matrix in text format (with values separated by a single space).

The default output is a unitig matrix whose values represent the fraction of the unitig's k-mers belonging to a sample. The `-r`/`--min-utg-frac` option allows to output a binary matrix which sets values to `1` when fraction is above the provided threshol and `0` otherwise.

Here is an example:

| Unitig ID | Sample 1 | Sample 2 | Sample 3 | Sample 4 | Sample 5 |
|-----------|----------|----------|----------|----------|----------|
| 0         | 0        | 1        | 0.23     | 0.3      | 1        |
| 1         | 1        | 1        |   0      | 0.8      | 0.4      |
| 2         | 0.47     | 0.2      |   1      |  1       | 0        |
| 3         | 0.8      | 1        | 0.78     |  1       | 0.81     |
| 4         | 0.79     | 1        |   1      | 0.87     | 0.89     |

With option `-r 0.8`, you would have:

| Unitig ID | Sample 1 | Sample 2 | Sample 3 | Sample 4 | Sample 5 |
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
