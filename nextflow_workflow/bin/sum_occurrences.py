#!/usr/bin/env python3
import argparse
import json
import sys
import typing
import numpy as np
from os import path
from typing import *


INDENT_VAL = 4


def overwrite_check(file_path: str, force: bool) -> None:
    if path.isfile(file_path) and not force:
        raise FileExistsError(
            f"Output file {file_path} already exists- either move or delete this file or enable --force")


def read_json(json_file: typing.TextIO, verbose: bool) -> Dict[str, np.ndarray]:
    occurrence_dict = {}

    # json.load() will return the .json as a Python dict
    json_dict = json.load(json_file)
    for vir_name, occurrence_list in json_dict.items():
        # casting the occurrence list to a numpy array will allow for fast element-wise addition later
        occurrence_dict[vir_name] = np.array(occurrence_list)

    return occurrence_dict


def read_json_from_path(input_path: str, verbose: bool) -> Dict[str, np.ndarray]:
    with open(input_path, "r") as json_file:
        return read_json(json_file, verbose)


def sum_occurrences(dict_list: List[Dict[str, np.ndarray]], verbose: bool) -> Dict[str, np.ndarray]:
    summed_dict = {}
    counter = 0

    # iterate through all individual occurrence dicts. For each dict, iterate through its keys and values, saving them
    # in a summed dictionary where the values for each virus are summed together
    for indv_dict in dict_list:
        for vir_name, occurrence_array in indv_dict.items():
            if vir_name in summed_dict:
                summed_dict[vir_name] += occurrence_array
            else:
                summed_dict[vir_name] = occurrence_array

        counter += 1

    if verbose:
        print(f"{counter} input individual occurrence count .jsons summed")

    return summed_dict


def write_occurrence_json(json_file: typing.TextIO, occurrence_dict: Dict[str, np.ndarray]) -> None:
    # convert numpy arrays back to lists so they'll work with the json package
    for key, array in occurrence_dict.items():
        occurrence_dict[key] = array.tolist()

    json_file.write(json.dumps(occurrence_dict, indent=INDENT_VAL))


def write_occurrence_json_from_path(json_path: str, occurrence_dict: Dict[str, np.ndarray], force: bool) -> None:
    overwrite_check(json_path, force)

    with open(json_path, "w") as occ_json:
        write_occurrence_json(occ_json, occurrence_dict)
        occ_json.close()


def parse_args(sys_args: list) -> argparse.Namespace:
    default_eval = 1e-5
    parser = argparse.ArgumentParser(sys_args, description="Accepts file with paths to VIBES integration .json files, then sums up occurrence counts for each reference with hits")
    parser.add_argument("output_path", type=str, help="Path to output .json file containing a summed nucleotide occurrence count for each reference virus that appears at least once in the input .jsons")
        # nargs="*" tells parsearg that the input type will be a list with at least 0 items
    parser.add_argument("json_path_list", type=str, nargs="*", help="Path to text file where each line is a path to a VIBES integration .json")
    parser.add_argument("--verbose", action="store_true", help="Print additional information useful for debugging")
    parser.add_argument("--force", action="store_true", help="If output file already exists, overwrite it")

    return parser.parse_args()

def _main():
    args = parse_args(sys.argv[1:])
    json_paths_list = args.json_path_list
    output_path = args.output_path
    verbose = args.verbose
    force = args.force

    dict_list = []

    for json_path in json_paths_list:
        json_path = json_path.strip("[],")
        individual_occurrences_dict = read_json_from_path(json_path, verbose)
        dict_list.append(individual_occurrences_dict)

    summed_occurrences_dict = sum_occurrences(dict_list, verbose)

    write_occurrence_json_from_path(output_path, summed_occurrences_dict, force)



    # The gist: read in every integration .json from the input file. Parse each file, building a dictionary of viruses
    # and occurrence counts. For each key in the dictionary, write a ,json containing the summed count.


if __name__ == "__main__":
    _main()