#!/usr/bin/env python3
import argparse
import sys
import re
from typing import TextIO
from typing import *
from os import path
from pathlib import Path
import json


DEFAULT_PREFIX = "contig_id"
INDENT_VAL = 4


def to_file(output_path: str, text: str, verbose: bool) -> None:
    if verbose:
        print(f"Writing reverted text to {output_path}...")

    with open(output_path, "w") as output_file:
        # get rid of escape character \ that changes original sequence name to something different
        text_no_escape = text.encode("utf-8").decode('unicode-escape')
        output_file.write(text_no_escape)


def revert_contig_ids(input_file: TextIO, id_mapping_dict: Dict[str, str], verbose: bool) -> str:
    file_content = input_file.read() # read entire input file in one go (so we have to use re.sub() less)

    for key in id_mapping_dict.keys(): # key here will be the replacement ID
        if verbose:
            print(f"Reverting {key} to {id_mapping_dict[key]}...")

        file_content = re.sub(key, id_mapping_dict[key], file_content)

    return file_content


def revert_contig_ids_from_path(input_path: str, id_mapping_dict: Dict[str, str], verbose: bool) -> str:
    if verbose:
        print(f"Opening {input_path}...")
    with open(input_path, "r") as input_file:
        return revert_contig_ids(input_file, id_mapping_dict, verbose)


def load_map_dict_from_json(json_path: str, verbose: bool) -> Dict[str, str]:
    map_dict = {}

    if verbose:
        print(f"Loading IDs from {json_path}...")

    with open(json_path, "r") as json_file:
         map_dict = json.load(json_file)

    return map_dict


def dict_to_json(json_path: str, map_dict: Dict[str, str]) -> None:
    with open(json_path, "w") as json_file:
        json_file.write(json.dumps(map_dict, indent=INDENT_VAL))


def tuple_list_to_fasta(fasta_path: str, tuple_list: List[Tuple[str,str,str]]) -> None:
    with open(fasta_path, "w") as output_fasta:
        for entry in tuple_list:
            header = entry[1]
            sequence = entry[2]
            fasta_entry = f">{header}\n{sequence}\n"

            output_fasta.write(fasta_entry)


def change_contig_ids(fasta_tuple_list: List[Tuple[str,str,str]], prefix: str, verbose: bool) -> Tuple[List[Tuple[str,str,str]], Dict[str,str]]:
    replaced_tuple_list = []
    id_map_dict = {}
    index = 1

    for entry in fasta_tuple_list:
        original_id = entry[0]
        original_header = entry[1]
        sequence = entry[2]

        # generate replacement values. Sequence shouldn't be affected by this
        replacement_id = f"{prefix}_{index}"
        replacement_header = re.sub(original_id, replacement_id, original_header)

        index += 1

        replaced_tuple_list.append((replacement_id, replacement_header, sequence))
        id_map_dict[replacement_id] = original_id

        if verbose:
            print(f"Replacing {original_id} with {replacement_id}...")

    # the values in the dict have added backslashes to escape any weird characters- we need to get rid of them so they
    # match the original contig ID exactly
    for key in id_map_dict.keys():
        id_map_dict[key] = re.sub("\\\\", "", id_map_dict[key])
        
    return replaced_tuple_list, id_map_dict


def parse_fasta(fasta_file: TextIO) -> List[Tuple[str, str, str]]:
    parsed_tuple_list = []

    fasta_list = re.split('\n>', fasta_file.read())

    for entry in fasta_list:
        # remove leading and trailing > character, if present
        entry = entry.strip(">")
        header, sequence = entry.split("\n", 1)
        # grab the 'name,' or fasta header line up to the first whitespace character
        seq_id = re.escape(header.split()[0])

        entry_tuple = (seq_id, header, sequence)
        parsed_tuple_list.append(entry_tuple)

    return parsed_tuple_list


def parse_fasta_from_path(fasta_path: str, verbose: bool) -> List[Tuple[str, str, str]]:
    with open(fasta_path, "r") as fasta_file:
        if verbose:
            print(f"Opening {fasta_path}...")

        return parse_fasta(fasta_file)


def parse_args(sys_args: list) -> argparse.Namespace:
    # default value for setting and use in help messages
    default_prefix = DEFAULT_PREFIX
    parser = argparse.ArgumentParser(sys_args, description="Either changes contig IDs in a .fasta file to replacement "
                                                           "IDs or reverts replacement IDs to originals in any text file")
    subparsers = parser.add_subparsers(dest='change_mode', help="Set program to either rename original IDs to "
                                                                "replacement IDs or revert replacements IDs back to "
                                                                "original IDs")
    # add parser for rename mode, which will rename contigs with a replacement ID and preserve the original ID in a JSON
    # we also add arguments unique to the rename parser
    change_parser = subparsers.add_parser('rename', help=f"Replace original contig ID with replacement ID based on a "
                                                         f"prefix (default {default_prefix}) and index as follows: "
                                                         f"{default_prefix}_1")
    change_parser.add_argument("input_fasta", type=str, help="Input .fasta format file with contigs to be renamed. By "
                                                             "default, this file will be overwritten to replace the "
                                                             "contig IDs (set --output_fasta to set a different output "
                                                             "file instead")
    change_parser.add_argument("output_map_json", type=str, help="Output .json file containing a converted Python "
                                                             "dictionary that maps new IDs to original IDs")
    change_parser.add_argument("--prefix", type=str, help=f"Sets prefix for replacement IDs (default {default_prefix})",
                               default=default_prefix)
    change_parser.add_argument("--output_fasta", type=str, help="Optional output .fasta file. Will be a copy of the "
                                                                "original input file with different contig IDs",
                               default="")
    # add parser for revert mode, which will scan through input file line by line, reverting contig ID to original
    # we also add arguments unique to the revert parser
    revert_parser = subparsers.add_parser('revert', help="Reverts any instance of a replacement contig ID with the "
                                                         "original contig ID")
    revert_parser.add_argument("input_file", type=str, help="Any file containing text that includes a replacement "
                                                            "contig ID. Each instance of a replacement contig ID will "
                                                            "be reverted to the original contig ID")
    revert_parser.add_argument("input_map_json", type=str, help=".json file produced by 'rename' mode. Contains a "
                                                                "Python dictionary mapping replacement IDs to original "
                                                                "IDs")
    revert_parser.add_argument("--output_file", type=str, help="Optional output file. Will be a copy of the original "
                                                               "input file with different contig IDs", default="")
# subparsers have independent lists of arguments. To set common arguments, loop over the subparsers
    for name, subp in subparsers.choices.items():
        subp.add_argument("--verbose", action="store_true", help="Set program to report whenever it modifies an ID")

    return parser.parse_args()

def _main():
    # TODO: Explanations
    args = parse_args(sys.argv[1:]) # get command line args. Right now we can only grab arguments both have in common
    change_mode = args.change_mode # sets whether we rename contig IDs or revert to original names
    verbose = args.verbose # reports when we change something in either mode

    if change_mode == 'rename':
        # now parse rename mode only args
        # contains contigs with IDs (anything between > and the first whitespace) to be renamed
        fasta_file = args.input_fasta
        # maps original IDs to replacement IDs
        output_json = args.output_map_json
        # set prefix for replacement IDs
        prefix = args.prefix
        # optional output file path. Empty and evaluates to False unless set by user
        output_fasta = args.output_fasta


        # unless set by user, this should be the same as the input .fasta
        if not output_fasta:
            output_fasta = fasta_file

        # break fasta file into tuples containing (ID, header, sequence) where ID is anything between > and the first
        # whitespace character and header is the full line after >
        entry_tuple_list = parse_fasta_from_path(fasta_file, verbose)

        # change IDs and get dictionary mapping replacement IDs (key) to original ID (value)
        modified_entry_tuple_list, id_map_dict = change_contig_ids(entry_tuple_list, prefix, verbose)

        # write modified .fasta entries, ID mapping structure to files
        tuple_list_to_fasta(output_fasta, modified_entry_tuple_list)
        dict_to_json(output_json, id_map_dict)

    elif change_mode == 'revert':
        input_path = args.input_file # input file with replaced contig IDs, to be reverted to original IDs
        input_json_path = args.input_map_json # input JSON containing dictionary mapping replaced ID (key) to original
                                                # ID (value)
        output_path = args.output_file # optional output file to save reverted IDs to. By default, we overwrite input_path

        if not output_path:
            output_path = input_path

        id_mapping_dict = load_map_dict_from_json(input_json_path, verbose) # load input JSON to get ID mapping dict

        file_text = revert_contig_ids_from_path(input_path, id_mapping_dict, verbose) # load input file to replace IDs

        to_file(output_path, file_text, verbose) # write modified file text to output file (default: input file)


if __name__ == "__main__":
    _main()