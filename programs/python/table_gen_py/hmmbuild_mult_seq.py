import re
import sys
import argparse
import subprocess
from io import TextIOWrapper
from typing import *


def generate_hmmbuild_cmd(temp_fasta_path: str, temp_hmm_path: str, dna: bool, rna: bool, amino: bool, cpu_count: int) -> str:
    cmd = "hmmbuild"

    if dna:
        cmd = f"{cmd} --dna"
    elif rna:
        cmd = f"{cmd} --rna"
    elif amino:
        cmd = f"{cmd} --amino"
    else:
        pass

    cmd = f"{cmd} --cpu {cpu_count} {temp_hmm_path} {temp_fasta_path}"

    return cmd


def generate_hmm(temp_fasta_list: List[str], dna: bool, rna: bool, amino: bool, cpu_count: int, verbose):

    for temp_fasta_path in temp_fasta_list:
        temp_hmm_path = re.sub(r"\.fasta", ".hmm", temp_fasta_path)
        cmd = generate_hmmbuild_cmd(temp_fasta_path, temp_fasta_path, dna, rna, amino, cpu_count)
        do_cmd(cmd, verbose)


def generate_temp_fastas(fasta_file: TextIOWrapper, temp_folder: str) -> List[str]:
    temp_fasta_list = []
    fasta_list = re.split('\n>', fasta_file.read())
    increment = 1

    for entry in fasta_list:
        # remove leading and trailing > character, if present
        entry = re.sub(r"^>|>$", "", entry, re.M)
        # first line is header, while the rest is the sequence body. uses capture groups before and after the first
        # newline to capture both
        header, sequence = re.match(r"(.+?)\n(.+)", entry, re.S).groups()

        temp_file_path = f"{temp_folder}{increment}.fasta"
        temp_fasta_list.append(temp_file_path)

        with open(temp_file_path, "w") as temp_file:
            temp_file.write(f">{header.strip()}\n{sequence.strip()}")
            temp_file.close()

        increment += 1

    return temp_fasta_list


def generate_temp_fastas_from_path(fasta_path: str, temp_folder: str) -> List[str]:
    with open(fasta_path) as fasta_file:
        temp_fasta_list = generate_temp_fastas(fasta_file, temp_folder)

    return temp_fasta_list


def do_cmd(cmd: str, verbose: bool):
    if verbose:
        print(f"Running command: {cmd}")
    subprocess.run(cmd.split())


def parse_args(sys_args: list) -> argparse.Namespace:
    parser = argparse.ArgumentParser(sys_args, description="Accepts input .fasta file and generates HMM for each entry. Automatically runs hmmpress on output .hmm file")
    parser.add_argument("input_fasta", type=str, help="Input .fasta format file containing dna/rna/amino acid sequences")
    parser.add_argument("output_hmm", type=str, help="Path to output .hmm file. Output.hmm will be accompanied by auxiliary 'pressed' files")
    parser.add_argument("--temp_folder", type=str, default=None, help="Path to folder where temporary .fasta files will be created. These are automatically deleted before the program ends")
    parser.add_argument("--verbose", help="Prints information about commands used, how many .fasta entries have been hmmbuilt", action="store_true")
    parser.add_argument("--cpu", type=int, default=1, help="How many threads hmmbuild will use (i > 0)")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--dna", default=False, help="Specifies that .fasta entries are DNA seq. Mutually exclusive with --rna, --amino", action="store_true")
    group.add_argument("--rna", default=False, help="Specifies that .fasta entries are RNA seq. Mutually exclusive with --dna, --amino", action="store_true")
    group.add_argument("--amino", default=False, help="Specifies that .fasta entries are amino seq. Mutually exclusive with --dna, --rna", action="store_true")

    return parser.parse_args()


def _main():
    args = parse_args(sys.argv[1:])
    fasta_path = args.input_fasta
    hmm_path = args.output_hmm
    temp_folder = args.temp_folder
    verbose = args.verbose
    cpu_count = args.cpu
    dna = args.dna
    rna = args.rna
    amino = args.amino

    if cpu_count < 0:
        raise ValueError("--cpu must be used with an argument greater than or equal to 0")

    if not temp_folder:
        temp_folder = re.sub(r"(.+)\/.+?\..+?$", "\g<1>/", fasta_path)

    temp_fasta_list = generate_temp_fastas_from_path(fasta_path, temp_folder)
    generate_hmm(temp_fasta_list, dna, rna, amino, cpu_count, verbose)


if __name__ == "__main__":
    _main()