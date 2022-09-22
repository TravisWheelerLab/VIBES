import subprocess
import argparse
from os import walk, path
import re
import sys
from typing import *


def dom_tbl_gen(prophage_name: str, input_dir: str, hmm_path:str, gen_code: int, output_dir: str, verbose: bool) -> None:

    output_file = "%s/%s.domtbl" % (output_dir, prophage_name)
    input_file = "%s/%s.fasta" % (input_dir, prophage_name)

    # if file exists already, and --force wasn't used, print a warning and exit
    if path.isfile(output_file) and not force:
        print(f"\nWarning: File {output_file} already exists. If you wish to overwrite files, use --force")
        exit()

    else:
        # run hmmscant to create .domtbl of prophage
        command = ["hmmscant", "-c", gen_code, "--domtblout", output_file, hmm_path, input_file]
        do_cmd(command, verbose)


# credit to Viktor Kerkez in: https://stackoverflow.com/questions/18160078/how-do-you-write-tests-for-the-argparse-portion-of-a-python-module
def parse_args(sys_args: list)-> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_dir", help="Path to directory of .fasta files. Expected format is [sequenceName].fasta")
    parser.add_argument("hmm_db", help="Path to HMM database to run against .fasta files")
    parser.add_argument("output_dir", help="Path to directory to store output .domtbl files")
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
