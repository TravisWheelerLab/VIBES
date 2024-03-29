import argparse
from os import path
import pandas as pd
import subprocess
from typing import *
import sys

DEFAULT_FLANK_LENGTH = 15000


def parse_args(sys_args: list) -> argparse.Namespace:
    parser = argparse.ArgumentParser(sys_args,
                                     description="Scans VIBES .tsv files for hits that match element_name and uses seqkit subseq to grab a region centered on the hit from the original bacterial genome")
    parser.add_argument("input_tsv", type=str,
                        help="Input VIBES .tsv file containing information on hits in bacterial genome")
    parser.add_argument("bacterial_genome", type=str,
                        help="Bacterial genome in .fasta format, which input_tsv is based on")
    parser.add_argument("element_name", type=str,
                        help="Name of element the grabbed region will be centered on. Element_name must exactly match the .tsv Name field")
    parser.add_argument("output_fasta", type=str,
                        help="Path to output .fasta file, which will contain sequence of regions around hits matching with name matching element_name")
    parser.add_argument("--flank_length", type=int, default=DEFAULT_FLANK_LENGTH,
                        help=f"Length of each flank to be grabbed from each end of the element (default {DEFAULT_FLANK_LENGTH}). The total length of the grabbed region will be ~twice the flank_length. If a flank would extend past either end of the genome, it will be cut short to match the end of the genome")
    parser.add_argument("--verbose",
                        help="Prints seqkit subseq commands run by region_grabber.py",
                        action="store_true")
    parser.add_argument("--force", help="If output file already exists, overwrite it", action="store_true")

    return parser.parse_args()


def read_tsv(tsv_path: str) -> pd.DataFrame:
    tsv_data = pd.read_csv(tsv_path, sep="\t")
    return tsv_data


def get_element_info(tsv_data: pd.DataFrame, element_name: str) -> pd.DataFrame:
    # grab all rows with name == element_name
    matching_rows = tsv_data.loc[tsv_data["Name"] == element_name]
    # reduce retained columns to relevant info: Name, Query Name, bacterial start, bacterial end, bacterial length, +/-
    matching_rows = matching_rows[["Hit Name", "Query Sequence Name", "Match Start on Query Seq",
                                   "Match End on Query Seq", "Query Genome Length", "Strand"]]

    return matching_rows


def get_subseq_list(element_data: pd.DataFrame, flank_size: int, genome_path: str, verbose: bool) -> List[str]:
    subseq_list = []
    for row in element_data.itertuples():
        origin_seq_name = row[2]
        match_start = row[3]
        match_end = row[4]
        genome_total_length = row[5]
        strand = row[6]

        start_pos = None
        end_pos = None
        # check strand: if -, then match start and end should be flipped as the hit is inverted and starts at match_end
        if strand == "+":
            start_pos = match_start - flank_size
            end_pos = match_end + flank_size
        elif row[6] == "-":
            start_pos = match_end - flank_size
            end_pos = match_start + flank_size
        else:
            raise ValueError(f"Unexpected strand type indicator {row[6]: expected + or -}")

        if start_pos < 1:
            start_pos = 1

        if end_pos > genome_total_length: # if greater than genome length
            end_pos = genome_total_length

        # seqkit subseq retrieves subsequence from genome with seq name origin_seq_name, coordinates start_pos-end_pos
        cmd = ["seqkit", "subseq", genome_path, "--chr", origin_seq_name, "-r", f"{start_pos}:{end_pos}"]
        subseq_list.append(do_cmd(cmd, verbose).stdout)

    return subseq_list


def write_to_output(output_path: str, subseq_list: List[str], force: bool) -> None:
    if path.isfile(output_path) and not force:
        raise FileExistsError(
            f"Output file {output_path} already exists- either move or delete this file or enable --force")
    else:
        with open(output_path, "w") as output_file:
            for entry in subseq_list:
                entry = entry + "\n"
                output_file.write(entry)


def index_cleanup(genome_path: str, verbose: bool) -> None:
    # this function deletes the index file created by seqkit
    index_path = genome_path + ".seqkit.fai"
    cmd = ["rm", index_path]
    do_cmd(cmd, verbose)


def do_cmd(cmd: List[str], verbose: bool) -> subprocess.CompletedProcess:
    # double check that all elements are strings
    cmd = [str(e) for e in cmd]
    if verbose:
        verbose_cmd = " ".join(cmd)
        print(f"Running command: {verbose_cmd}")

    return subprocess.run(cmd, capture_output=True, text=True)


def _main():
    args = parse_args(sys.argv[1:])
    tsv_path = args.input_tsv
    genome_path = args.bacterial_genome
    element_name = args.element_name
    output_path = args.output_fasta
    flank_len = args.flank_length
    verbose = args.verbose
    force = args.force

    tsv_data = read_tsv(tsv_path)
    element_data = get_element_info(tsv_data, element_name)
    subseqs = get_subseq_list(element_data, flank_len, genome_path, verbose)

    if len(subseqs) > 0:
        write_to_output(output_path, subseqs, force)

    index_cleanup(genome_path, verbose)


if __name__ == "__main__":
    _main()
