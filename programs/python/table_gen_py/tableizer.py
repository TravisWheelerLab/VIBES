import argparse
import sys
import subprocess
from typing import *
from os import path
from os import remove


def generate_scanned_dfam(hmm_path: str, genome_path: str, dfam_dir: str, cpu_count: int, verbose: bool, force: bool):
    input_name = path.basename(hmm_path)
    dfam_path = f"{dfam_dir}/{input_name}.dfam"
    unscanned_dfam_path = f"{dfam_dir}/{input_name}_unscanned.dfam"

    # if user enabled --force, remove old output files before continuing. Always remove unscanned file if it exists
    if path.exists(unscanned_dfam_path):
        remove(unscanned_dfam_path)

    if path.exists(dfam_path):
        if force:
            remove(dfam_path)
        else:
            raise FileExistsError(f"Output file {dfam_path} already exists- either move or delete this file or enable --force")

    nhmmscan_cmd = ["nhmmscan", "--cpu", cpu_count, "--dfamtblout", unscanned_dfam_path, hmm_path, genome_path]
    do_cmd(nhmmscan_cmd, verbose)

    dfamscan_cmd = ["./dfamscan.pl", "--dfam_infile", unscanned_dfam_path, "--dfam_outfile", dfam_path]
    do_cmd(dfamscan_cmd, verbose)

    remove(unscanned_dfam_path)


def do_cmd(cmd: List[str], verbose: bool):
    # double check that all elements are strings
    for i, element in enumerate(cmd):
        cmd[i] = str(element)
    if verbose:
        verbose_cmd = " ".join(cmd)
        print(f"Running command: {verbose_cmd}")

    subprocess.run(cmd)


def parse_args(sys_args: list) -> argparse.Namespace:
    parser = argparse.ArgumentParser(sys_args, description="Accepts path to a HMM and a genome. Uses nhmmscan to produce output .dfam files, which are automatically"
                                                           " scanned by dfamscan.pl to resolve overlapping hits")
    parser.add_argument("hmm_path", type=str, help="Path to input .hmm file. The input .hmm file's directory must also contain auxiliary .hmm.h3f, "
                                              ".hmm.h3i, .hmm.h3m, and .hmm.h3p files, which are generated by hmmpress (automatically handled by hmmbuild_mult_seq.py")
    parser.add_argument("genome_path", type=str, help="Path to input genome in .fasta format")
    parser.add_argument("output_table_dir", type=str, help="Path to output directory of .dfam table files")
    parser.add_argument("--cpu", type=int, default=1, help="How many threads nhmmscan will use (i > 0)")
    parser.add_argument("--verbose", action="store_true", help="Print additional information useful for debugging")
    parser.add_argument("--force", action="store_true", help="If output file already exists, overwrite it")

    return parser.parse_args()


def _main():
    args = parse_args(sys.argv[1:])
    hmm_path = args.hmm_path
    genome_path = args.genome_path
    dfam_dir = args.output_table_dir
    cpu_count = args.cpu
    verbose = args.verbose
    force = args.force

    # check that inputs are legal
    if cpu_count < 0:
        raise ValueError("--cpu must be used with an argument greater than or equal to 0")

    if not path.isfile(hmm_path):
        raise FileNotFoundError(f"No such file: {hmm_path}")

    if not path.isfile(genome_path):
        raise FileNotFoundError(f"No such file: {genome_path}")

    if not path.isdir(dfam_dir):
        raise NotADirectoryError(f"No such folder: {dfam_dir}")

    generate_scanned_dfam(hmm_path, genome_path, dfam_dir, cpu_count, verbose, force)


if __name__ == "__main__":
    _main()