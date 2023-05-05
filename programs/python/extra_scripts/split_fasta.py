import argparse
import sys
from typing import Dict, TextIO
import re


def generate_individual_fastas(fasta_file: TextIO, output_dir: str):
    temp_fasta_dict = {}
    fasta_list = re.split('\n>', fasta_file.read())
    # TODO: explain why not >

    for entry in fasta_list:
        # remove leading and trailing > character, if present
        entry = entry.strip(">")
        header, sequence = entry.split("\n", 1)
        # TODO: explain magic number
        fasta_name = header.split()[0]
        # grab the 'name,' or fasta header line up to the first whitespace character
        seq_name = re.escape(header.split()[0])

        temp_file_path = f"{output_dir}/{fasta_name}.fasta"
        temp_fasta_dict[temp_file_path] = seq_name
        with open(temp_file_path, "w") as temp_file:
            temp_file.write(f">{header.strip()}\n{sequence.strip()}")
            temp_file.close()



def parse_args(sys_args: list) -> argparse.Namespace:
    parser = argparse.ArgumentParser(sys_args, description="Break a multi-sequence .fasta file into multiple files "
                                                           "containing once sequence each. Each file will be named "
                                                           "after the first word in its header line (whatever is "
                                                           "between the '>' character and the first whitespace "
                                                           "character")
    parser.add_argument("input_fasta", type=str, help="Path to input .fasta format file containing multiple entries "
                                                      "or sequences separated by a '>' character")
    parser.add_argument("fasta_dir", type=str, help="Path to parent output directory, where individual sequence files "
                                                    "will live")

    return parser.parse_args()


def _main():
    # get arguments from argparse
    args = parse_args(sys.argv[1:])
    input_fasta = args.input_fasta
    output_dir = args.fasta_dir

    with open(input_fasta) as fasta_file:
        generate_individual_fastas(fasta_file, output_dir)


if __name__ == "__main__":
    _main()