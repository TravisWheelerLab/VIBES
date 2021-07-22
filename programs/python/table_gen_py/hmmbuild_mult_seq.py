import re
import sys
import argparse
import subprocess
from typing import TextIO
from typing import *
from os import path
from os import remove


def hmmpress_output(output_path, verbose):
    do_cmd(["hmmpress", output_path], verbose)


def combine_hmms(temp_hmm_list: List[str], output_hmm_path: str, force: bool):
    if path.exists(output_hmm_path):
        if force:
            # TODO: this doesn't clean up hmmpress aux files
            remove(f"{output_hmm_path}")
        else:
            raise FileExistsError(f"Output file {output_hmm_path} already exists- either move or delete this file or enable --force")

    with open(output_hmm_path, "a") as output_file:
        for temp_hmm_path in temp_hmm_list:
            with open(temp_hmm_path) as temp_file:
                temp_contests = temp_file.read()
                output_file.write(f"{temp_contests}\n")


def generate_hmmbuild_cmd(temp_fasta_path: str, temp_hmm_path: str, dna: bool, rna: bool, amino: bool, cpu_count: int, seq_name: str) -> str:
    seq_type = ""

    if dna:
        seq_type = "--dna"
    elif rna:
        seq_type = "--rna"
    elif amino:
        seq_type = "--amino"
    else:
        pass

    cmd = ["hmmbuild", "--cpu", cpu_count, "-n",  seq_name, seq_type, temp_hmm_path, temp_fasta_path]

    return cmd


def generate_hmm(temp_fasta_dict: Dict[str, str], dna: bool, rna: bool, amino: bool, cpu_count: int, verbose) -> List[str]:
    temp_hmm_list = []

    for temp_fasta_path, seq_name in temp_fasta_dict.items():
        temp_hmm_path = f"{path.splitext(temp_fasta_path)[0]}.hmm"
        temp_hmm_list.append(temp_hmm_path)
        cmd = generate_hmmbuild_cmd(temp_fasta_path, temp_hmm_path, dna, rna, amino, cpu_count, seq_name)
        do_cmd(cmd, verbose)

    return temp_hmm_list


def generate_temp_fastas(fasta_file: TextIO, temp_folder: str) -> Dict[str, str]:
    temp_fasta_dict = {}
    fasta_list = re.split('\n>', fasta_file.read())
    increment = 1

    for entry in fasta_list:
        # remove leading and trailing > character, if present
        entry = entry.strip(">")
        # first line is header, while the rest is the sequence body. uses capture groups before and after the first
        # newline to capture both
        header, sequence = entry.split("\n", 1)
        # grab the 'name,' or fasta header line up to the first whitespace character
        seq_name = re.escape(header.split()[0])

        temp_file_path = f"{temp_folder}temp{increment}.fasta"
        temp_fasta_dict[temp_file_path] = seq_name
        with open(temp_file_path, "w") as temp_file:
            temp_file.write(f">{header.strip()}\n{sequence.strip()}")
            temp_file.close()

        increment += 1

    return temp_fasta_dict


def generate_temp_fastas_from_path(fasta_path: str, temp_folder: str) -> Dict[str, str]:
    with open(fasta_path) as fasta_file:
        temp_fasta_dict = generate_temp_fastas(fasta_file, temp_folder)

    return temp_fasta_dict


def do_cmd(cmd: List[str], verbose: bool):
    # double check that all elements are strings
    for i, element in enumerate(cmd):
        cmd[i] = str(element)
    print(cmd)
    if verbose:
        verbose_cmd = " ".join(cmd)
        print(f"Running command: {verbose_cmd}")

    subprocess.run(cmd)


def parse_args(sys_args: list) -> argparse.Namespace:
    parser = argparse.ArgumentParser(sys_args, description="Accepts input .fasta file and generates HMM for each entry. Automatically runs hmmpress on output .hmm file")
    parser.add_argument("input_fasta", type=str, help="Input .fasta format file containing dna/rna/amino acid sequences")
    parser.add_argument("output_hmm", type=str, help="Path to output .hmm file. Output.hmm will be accompanied by auxiliary 'pressed' files")
    parser.add_argument("--temp_folder", type=str, default=None, help="Path to folder where temporary .fasta files will be created. These are automatically deleted before the program ends."
                                                                      "If no folder is specified, temporary files are stored in the directory that the output file will live in")
    parser.add_argument("--verbose", help="Prints information about commands used, how many .fasta entries have been hmmbuilt", action="store_true")
    parser.add_argument("--force", help="If output file already exists, overwrite it", action="store_true")
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
    force = args.force
    cpu_count = args.cpu
    dna = args.dna
    rna = args.rna
    amino = args.amino

    if cpu_count < 0:
        raise ValueError("--cpu must be used with an argument greater than or equal to 0")

    if not temp_folder:
        # this regex statement should grab the path of the input file up to its last / character (the path to the input
        # .fasta file's directory). we use this as the temporary file folder unless an alternative has been provided by
        # the user
        temp_folder = re.sub(r"(.+)\/.+?\..+?$", "\g<1>/", fasta_path)

    temp_fasta_dict = generate_temp_fastas_from_path(fasta_path, temp_folder)
    temp_hmm_list = generate_hmm(temp_fasta_dict, dna, rna, amino, cpu_count, verbose)

    combine_hmms(temp_hmm_list, hmm_path, force)
    hmmpress_output(hmm_path, verbose)

    for file_path in temp_fasta_dict.keys():
        remove(file_path)

    for file_path in temp_hmm_list:
        remove(file_path)

if __name__ == "__main__":
    _main()