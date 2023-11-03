#!/usr/bin/env python3
import subprocess
import argparse
from os import walk, path
import re
import sys
from typing import *


def dom_tbl_gen(prophage_name: str, input_path: str, hmm_db_path:str, gen_code: int, output_path: str, verbose: bool) -> None:
    # if file exists already, and --force wasn't used, print a warning and exit
    if path.isfile(output_path) and not force:
        print(f"\nWarning: File {output_path} already exists. If you want to overwrite files, use --force")
        exit()

    else:
        # run hmmscant to create .domtbl of prophage
        # TODO: Replace with FraHMMER
        command = ["frahmmer", "-c", gen_code, "--tblout", output_path, input_path, hmm_db_path]
        do_cmd(command, verbose)


def parse_args(sys_args: list)-> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_dir", help="Path to directory of .fasta files. Expected format is [sequenceName].fasta")
    parser.add_argument("hmm_db", help="Path to HMM database to run against .fasta files")
    parser.add_argument("output_path", help="Path to directory to store output .domtbl files")
    parser.add_argument("--genetic_code", type=int, default=1, help="Instruct hmmscant to use alt genetic code of NCBI translation table. Default is 1 (standard)")
    parser.add_argument("--force", help="If output files already exist, overwrite them", action="store_true")
    parser.add_argument("--verbose", help="Print commands run by generate_domtbls.py", action="store_true")
    return parser.parse_args()


def do_cmd(cmd: List[str], verbose: bool) -> None:
    # double check that all elements are strings
    cmd = [str(e) for e in cmd]
    if verbose:
        verbose_cmd = " ".join(cmd)
        print(f"Running command: {verbose_cmd}")

    subprocess.run(cmd)


if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    input_dir = args.fasta_dir
    hmm_path = args.hmm_db
    output_dir = args.output_dir
    gen_code = args.genetic_code
    force = args.force
    verbose = args.verbose

    # credit to pycruft in https://stackoverflow.com/questions/3207219/how-do-i-list-all-files-of-a-directory for code for grabbing file paths
    file_paths = []
    for (dir_path, dir_names, files) in walk(input_dir):
        file_paths.extend(files)

    file_paths.sort()

    for file_path in file_paths:
        regMatch = re.match(r'(.+?)\.fasta', file_path)
        # extract name from file path
        prophage_name = regMatch.group(1)

        dom_tbl_gen(prophage_name, input_dir, hmm_path, gen_code, output_dir, verbose)
