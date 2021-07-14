import re
import sys
import argparse
import subprocess
from io import TextIOWrapper


def open_input_fasta(file_path:str):
    with open(file_path, 'r') as fasta_file:
        generate_temp_fasta(fasta_file)


def generate_temp_fasta(fasta_file: TextIOWrapper):
    fasta_list = re.split('\n>', fasta_file.read())
    print(fasta_list[1])



def do_cmd(cmd: str, verbose: bool):
    if verbose:
        print(f"Running command: {cmd}")
    subprocess.run(cmd.split())


def parse_args(sys_args: list) -> argparse.Namespace:
    parser = argparse.ArgumentParser(sys_args, description="Accepts input .fasta file and generates HMM for each entry. Automatically runs hmmpress on output .hmm file")
    parser.add_argument("input_fasta", type=str, help="Input .fasta format file containing dna/rna/amino acid sequences")
    parser.add_argument("output_hmm", type=str, help="Path to output .hmm file. Output.hmm will be accompanied by auxiliary 'pressed' files")
    parser.add_argument("--verbose", help="Prints information about commands used, how many .fasta entries have been hmmbuilt", action="store_true")
    parser.add_argument("--cpu", type=int, default=1, help="How many threads hmmbuild will use (i > 0)")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--dna", default=False, help="Specifies that .fasta entries are DNA seq", action="store_true")
    group.add_argument("--rna", default=False, help="Specifies that .fasta entries are RNA seq", action="store_true")
    group.add_argument("--amino", default=False, help="Specifies that .fasta entries are amino seq", action="store_true")

    return parser.parse_args()


def _main():
    args = parse_args(sys.argv[1:])
    fasta_path = args.input_fasta
    hmm_path = args.output_hmm
    verbose = args.verbose
    cpu_count = args.cpu
    dna_seq = args.dna
    rna_seq = args.rna
    amino_seq = args.amino

    if cpu_count < 0:
        raise ValueError("--cpu must be used with an argument greater than or equal to 0")

    open_input_fasta(fasta_path)


if __name__ == "__main__":
    _main()