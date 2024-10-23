#!/usr/bin/env python3

import os, sys, shutil, argparse, logging
import pathlib

logger = logging.getLogger()

MUSET_VERSION="0.4.1"

MUSET_BIN_PATH = pathlib.Path(__file__).resolve()
MUSET_BIN_DIR = MUSET_BIN_PATH.parent


def init_logging(debug:bool) -> None:
    global logger
    log_formatter = logging.Formatter('[{asctime}] {message}', datefmt='%Y-%m-%d %H:%M:%S', style='{')
    stream_handler = logging.StreamHandler(stream=sys.stderr)
    stream_handler.setFormatter(log_formatter)
    logger.setLevel(logging.DEBUG if debug else logging.INFO)
    logger.addHandler(stream_handler)

    # h_stdout = logging.StreamHandler(stream=sys.stdout)
    # h_stderr = logging.StreamHandler(stream=sys.stderr)
    # h_stderr.addFilter(lambda record: record.levelno >= logging.WARNING)
    # logging.basicConfig(
    #     level = args.level.upper(),
    #     format = "%(asctime)s [%(levelname)s] %(message)s",
    #     handlers = [
    #         logging.FileHandler(args2logname(args)),
    #         h_stdout,
    #         h_stderr
    #     ]
    # )


def main() -> int:

    def check_int_range(value, min, max):
        value = int(value)
        if value < min or value > max:
            raise argparse.ArgumentTypeError(f'value should be in the range [{min}, {max}].')
        return value

    def check_float_range(value, min, max):
        value = float(value)
        if value < min or value > max:
            raise argparse.ArgumentTypeError(f'value should be in the range [{min}, {max}].')
        return value

    parser = argparse.ArgumentParser(prog='muset', 
                                     description='A pipeline for building an abundance unitig matrix from a list of FASTA/FASTQ files',
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     add_help=False)
    parser.add_argument('fof', metavar='INPUT_FILE', nargs='?', help='kmtricks-like input file (see README.md). It is ignored if -i option is used.')
    
    # MAIN OPTIONS
    main_args = parser.add_argument_group('main options')
    main_args.add_argument('-i', dest='in_matrix', metavar='PATH', required=False, help='run the pipeline with a previosuly computed matrix.')
    main_args.add_argument('-o', dest='out_dir', metavar='PATH', required=False, type=str, default='output', help='output directory. {%(default)s}')
    main_args.add_argument('-k', dest='kmer_size', metavar='INT', required=False, type=lambda val:check_int_range(val,8,63), default=31, help='k-mer size. [8, 63] {%(default)s}')
    main_args.add_argument('-m', dest='mini_size', metavar='INT', required=False, type=lambda val:check_int_range(val,4,15), default=15, help='minimizer size. [4, 15] {%(default)s}')
    main_args.add_argument('-a', dest='min_abundance', metavar='INT', required=False, type=int, default=2, help='min abundance to keep a k-mer. {%(default)s}')
    main_args.add_argument('-l', dest='min_utg_len', metavar='INT', required=False, type=int, help='min length to keep a unitig in the output matrix. {2k-1}')
    main_args.add_argument('-s', dest='write_utg_seq', action='store_true', help='write the unitig sequence (instead of the identfier) in the first column of the output matrix.')

    # ADVANCED OPTIONS
    filter_args = parser.add_argument_group('filtering options')
    filter_args.add_argument('--no-kmer-filter', dest='disable_filter', action='store_true', help='disable filtering of k-mer matrix rows before unitig construction.')
    filter_args.add_argument('-f', dest='min_frac_absent', metavar='FLOAT', required=False, type=lambda val:check_float_range(val,0,1), default=0.1, help='fraction of samples from which a k-mer should be absent. [0,1] {%(default)s}')
    filter_args.add_argument('-F', dest='min_frac_present', metavar='FLOAT', required=False, type=lambda val:check_float_range(val,0,1), default=0.1, help='fraction of samples in which a k-mer should be present. [0,1] {%(default)s}')
    filter_args.add_argument('-n', dest='min_nb_absent', metavar='INT', required=False, type=int, help='minimum number of samples from which a k-mer should be absent (overrides -f).')
    filter_args.add_argument('-N', dest='min_nb_present', metavar='INT', required=False, type=int, help='minimum number of samples in which a k-mer should be present (overrides -F).')
    
    # OTHER OPTIONS
    other_args = parser.add_argument_group('other options')
    other_args.add_argument('-t','--threads', dest='nb_threads', metavar='INT', type=int, default=4, help='number of threads. {%(default)s}')
    other_args.add_argument('-h','--help', action='help', help='show this help message and exit')
    other_args.add_argument('-v','--version', action='version', version=f'muset {MUSET_VERSION}', help='show version number and exit')
    other_args.add_argument('--verbose', dest='debug', action='store_true', help='verbose output')
    opt = parser.parse_args()

    init_logging(opt.debug)

    # check companion executables exist

    kmat_tools = shutil.which("kmat_tools")
    if not kmat_tools:
        logger.error(f"kmat_tools binary was not found in path: {MUSET_BIN_DIR}")

    muset_kmtricks = shutil.which("muset-kmtricks")
    if not muset_kmtricks:
        logger.error(f"muset-kmtricks binary was not found in path: {MUSET_BIN_DIR}")

    # check input files

    if opt.in_matrix and (not os.path.isfile(opt.in_matrix) or os.path.getsize(opt.in_matrix)):
        logger.error(f"Input matrix \"{opt.in_matrix}\" is not a file or is empty")
        return 1

    if not opt.in_matrix and not opt.fof:
        logger.error(f"Input file missing.\nTry 'muset -h' for more information.")
        return 1
    
    if not opt.in_matrix and not os.path.isfile(opt.fof):
        logger.error(f"Input file \"{opt.fof}\" does not exist.")
        return 1
    
    # check minimizer size is smaller than k-mer size
    if opt.mini_size >= opt.kmer_size:
        logger.error(f"Minimizer size (-m) value ({opt.mini_size}) must be smaller than k-mer size (-k) value ({opt.kmer_size}).")
        return 1
    
    # Check if -n/-N or -f/-F have been properly provided or set default filter
    if (opt.min_nb_absent and not opt.min_nb_present) or (not opt.min_nb_absent and opt.min_nb_present):
        logger.error(f"Minimizer size (-m) value ({opt.mini_size}) must be smaller than k-mer size (-k) value ({opt.kmer_size}).")

    # set unitig length to 2k-1 if no custom value is provided
    if not opt.min_utg_len:
        opt.min_utg_len = 2 * opt.kmer_size - 1

    logger.info(f"muset v{MUSET_VERSION}")
    logger.info(f"k-mer size (-k): {opt.kmer_size}")
    logger.info(f"Number of threads (-t): {opt.nb_threads}")
    logger.info(f"Minimum unitig size (-l): {opt.min_utg_len}")
    logger.info(f"Minimizer length (-m): {opt.mini_size}")
    logger.info(f"Minimum abundance (-a): {opt.min_abundance}")
    logger.info(f"Skip matrix construction (-i): {opt.in_matrix is not None}")
    logger.info(f"Input file: {opt.fof}")
    logger.info(f"Output directory (-o): {opt.out_dir}")
    logger.info(f"Write unitig sequence instead of identifier (-s): {opt.write_utg_seq}")
    
    logger.info(f"Filter settings:")

    if opt.min_nb_absent:
        logger.info(f"Min number of samples absent: {opt.min_nb_absent}")
    else:
        logger.info(f"Min fraction of samples absent: {opt.min_frac_absent}")
    
    if opt.min_nb_present:
        logger.info(f"Min number of samples present: {opt.min_nb_present}")
    else:
        logger.info(f"Min fraction of samples present: {opt.min_frac_present}")


    return 0


if __name__ == "__main__":
    os.environ["PATH"] = str(MUSET_BIN_DIR) + os.pathsep + os.environ["PATH"]
    sys.exit(main())