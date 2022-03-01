import argparse
import sys
import subprocess
from typing import *
from os import path
from os import remove


def generate_scanned_dfam(hmm_path: str, genome_path: str, output_path: str, dfamscan_path: str, cpu_count: int, verbose: bool, force: bool):
    unscanned_dfam_path = path.join(path.dirname(output_path), f"{path.basename(output_path)}_unscanned.dfam")

    if path.exists(unscanned_dfam_path):
        remove(unscanned_dfam_path)

    if path.exists(output_path):
        if force:
            remove(output_path)
        else:
            raise FileExistsError(f"Output file {output_path} already exists- either move or delete this file or enable --force")

    nhmmscan_cmd = ["nhmmscan", "--cpu", cpu_count, "--dfamtblout", unscanned_dfam_path, hmm_path, genome_path]
    do_cmd(nhmmscan_cmd, verbose)

    dfamscan_cmd = [dfamscan_path, "--dfam_infile", unscanned_dfam_path, "--dfam_outfile", output_path]
    do_cmd(dfamscan_cmd, verbose)

    remove(unscanned_dfam_path)


def do_cmd(cmd: List[str], verbose: bool):
    # double check that all elements are strings
    cmd = [str(e) for e in cmd]
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
    parser.add_argument("output_table", type=str, help="Path to output, scanned .dfam file")
    parser.add_argument("--dfamscan_path", type = str, default="./dfamscan.pl", help="Path to dfamscan.pl library script. Only needs to be set if you're not running from same dir as this script")
    parser.add_argument("--cpu", type=int, default=1, help="How many threads nhmmscan will use (i > 0)")
    parser.add_argument("--verbose", action="store_true", help="Print additional information useful for debugging")
    parser.add_argument("--force", action="store_true", help="If output file already exists, overwrite it")

    return parser.parse_args()


def _main():
    args = parse_args(sys.argv[1:])
    hmm_path = args.hmm_path
    genome_path = args.genome_path
    output_path = args.output_table
    dfamscan_path = args.dfamscan_path
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

    #if not path.isfile(path.dirname(dfam_path)):
    #    raise NotADirectoryError(f"No such folder: {path.dirname(dfam_path)}")

    generate_scanned_dfam(hmm_path, genome_path, output_path, dfamscan_path, cpu_count, verbose, force)


if __name__ == "__main__":
    _main()