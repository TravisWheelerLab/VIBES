# gene_count_filter.py
# Description: This script ingests output VIBES TSV files (integration 
# identification and viral gene annotation) and counts how many genes occur in
# each match, filtering output lines that fail to meet a user-specified
# threshold. Based on a seperate, percentage-based threshold, we double check 
# that a meaningful chunk of each counted gene actually falls in a match.
# This program is intended to be used after the VIBES pipeline has already 
# completed execution.
# Author: Conner Copeland
# Contact: connercopeland01@gmail.com
# Created on: 2024-01-24

import argparse
import matplotlib.pyplot as plt
import os
from pathlib import Path
import sys
from typing import List, Dict, Tuple

# Minimum number of genes a match must contain to pass through the filter
DEFAULT_GENE_COUNT_THRESHOLD = 2

# Minimum percentage of a gene covered by a match for it to be counted towards passing the filter
DEFAULT_GENE_COVERAGE_THRESHOLD = .2


def plot_histogram(histogram_dict: Dict[int, int]) -> None:
    """
    Uses matplotlib to generate a histogram with genes per match information.

    Args:
         -histogram_dict (Dict[int, int]): Dictionary containing keys that represent a number of genes in a match
            and values that represent how many matches contained that many genes.

    Returns:
        None
    """
    plt.bar(histogram_dict.keys(), histogram_dict.values())
    plt.xlabel("Genes per Match")
    plt.ylabel("Matches")
    plt.title("Distribution of Genes per Match")


def write_lines_to_output(output_path: str, output_content: str, force: bool) -> None:
    """
    Opens output_path as a file and writes output_content to it. If the file exists, --force must be enabled
    """
    # If this file exists, unless force is enabled, don't write to the file
    if not Path(output_path).exists() or force:
        with open(output_path, "w") as output_file:
            output_file.write(output_content)
    else:
        raise FileExistsError(
            f"Output file {output_path} already exists- either move or delete this file or enable --force")


def sufficient_overlap(match_st: int, match_en: int, gene_st: int, gene_en: int, minimum_coverage: float) -> bool:
    """
    Returns True if the gene overlaps enough with the match to pass a threshold determined by
        DEFAULT_GENE_COVERAGE_THRESHOLD; otherwise, returns False.

    Args:
        - match_st (int): Start position of the match.
        - match_en (int): End position of the match.
        - gene_st (int): Start position of the gene.
        - gene_en (int): End position of the gene.
        -minimum_coverage (float): float representation of the minumum percentage by which the match and gene must
            overlap to be reported

    Returns:
        - sufficiently_overlaps (bool): True when overlap passes threshold, False otherwise.
    """
    # we've ordered the Tuple such that the smaller value is always gene_st, and the larger is always gene_en
    coverage_factor = (min(match_en, gene_en) - max(match_st, gene_st) / (gene_en - gene_st))

    return coverage_factor >= minimum_coverage


def gene_in_match(match_st: int, match_en: int, gene_st: int, gene_en: int) -> bool:
    """
    Returns True if the gene occurs at all within the bounds of the match; otherwise, returns False.

    Args:
        - match_st (int): Start position of the match.
        - match_en (int): End position of the match.
        - gene_st (int): Start position of the gene.
        - gene_en (int): End position of the gene.

    Returns:
        overlaps (bool): True when the match and gene overlap at all, False otherwise.
    """
    # If the end of the gene occurs after the start of the match, and the start of the gene occurs before the end of the
    # match, there must be at least some overlap.
    return (gene_en > match_st) and (gene_st < match_en)


def report_genes_in_match(match_tuple: Tuple[str, int, int], gene_dict: Dict[str, List[Tuple[int,int]]],
                          minimum_coverage: float) -> int:
    """
    Accepts a tuple containing match information (from one match on one bacteria genome) and viral gene annotation
    information (all query viruses). Detects genes overlapping with a match and returns how many genes do so.

    Args:
        - match_tuple (Tuple[str, int, int]): Contains information on one match: (match_name, query_st, query_en)
        - gene_list (Dict[str, List[Tuple[str, str]]]): Contains gene position information for each query virus:
            Dict[query_name, List[Tuple[gene_st, gene_en]]
        -minimum_coverage (float): float representation of the minumum percentage by which the match and gene must
            overlap to be reported

    Returns:
        - gene_in_match (int): Number of genes that overlap with match
    """
    genes_in_match = 0
    # use query virus name to identify correct gene position list
    gene_list = gene_dict[match_tuple[0]]

    match_st = match_tuple[1]
    match_en = match_tuple[2]

    for gene_tuple in gene_list:
        gene_st = gene_tuple[0]
        gene_en = gene_tuple[1]

        if gene_in_match(match_st, match_en, gene_st, gene_en) and sufficient_overlap(match_st, match_en, gene_st,
                                                                                      gene_en, minimum_coverage):
            genes_in_match += 1

    return genes_in_match


def generate_match_tuple(line: str) -> Tuple[str, int, int]:
    """
    Accepts a line from an output VIBES TSV file containing viral match information. Returns tuple containing
    the nearest query name and start and end coordinate of a match to a query virus on a bacterial genome.

    Args:
    - line (str): Line describing match in VIBES output bacterial_integration TSV

    Returns:
    - match_tuple (Tuple[str, int,int]): Tuple containing query name and start and end coordinates for
     a match.
    """
    match_tuple = ()
    split_line = line.split()

    # grab query name from TSV
    query_name = split_line[0]
    # grab gene's target (virus) start coord from TSV
    start_coord = int(split_line[5])
    # grab gene's target (virus) end coord from TSV
    end_coord = int(split_line[6])

    # if the end coordinate < start coordinate, the gene in on the negative strand. We don't care about that
    # here, so we just flip the coordinates- detecting overlap between matches and genes will work fine
    if start_coord < end_coord:
        match_tuple = (query_name, start_coord, end_coord)
    else:
        match_tuple = (query_name, end_coord, start_coord)

    return match_tuple


def build_query_virus_gene_dict(gene_annotation_paths: List[str]) -> Dict[str, List[Tuple[int, int]]]:
    """
    Builds a Dictionary where each key is the name of a query virus and values are a List of Tuples, each Tuple
    containing the start and end coordinate of an annotated gene on the associated query viral genome.

    Args:
    - gene_annotation_paths (List[str]): List of paths to viral gene annotation TSV files generated by VIBES

    Returns:
    - gene_dict (Dict[str, List[Tuple[int,int]]]): Dictionary where each key is the name of a query virus and
    values are a List of Tuples, each Tuple containing the start and end coordinate of an annotated gene on the
    associated query viral genome
    """
    gene_dict = {}

    for annotation_path in gene_annotation_paths:
        gene_list = generate_gene_list(annotation_path)
        # annotation TSVs are named after their viruses. With Path().stem, we get the file name without preceding
        # path or extensions
        gene_dict[Path(annotation_path).stem] = gene_list

    return gene_dict



def list_files_in_dir(dir_path: str) -> List[str]:
    """
    Returns a list containing all files in a directory.

    Args:
    - dir_path (str): Path to the directory containing files that we want.

    Returns:
    - file_list (List[str]): List of paths to files in the input directory.
    """
    file_list = []

    file_iterator = os.scandir(dir_path)

    for path in file_iterator:
        if path.is_file():
            file_list.append(path.path)


    return file_list

def generate_gene_list(annotation_tsv_path: str) -> List[Tuple[int, int]]:
    """
    Accepts path to TSV file containing viral gene annotation information. Returns a list of tuples, each containing
    the start and end coordinate of a gene.

    Args:
    - annotation_tsv_path (str): Path pointing to a VIBES output viral gene annotation TSV

    Returns:
    - gene_list (List[Tuple[int,int]]): List of tuples containing start and end coordinates for genes. Tuples occur in
    order based on start coordinate.
    """
    gene_list = []
    with open(annotation_tsv_path) as annotation_tsv_file:
        for line in annotation_tsv_file:
            if line[0] == '#':
                pass
            else:
                split_line = line.split()

                # grab gene's target (virus) start coord from TSV
                start_coord = int(split_line[10])
                # grab gene's target (virus) end coord from TSV
                end_coord = int(split_line[11])

                # if the end coordinate < start coordinate, the gene in on the negative strand. We don't care about that
                # here, so we just flip the coordinates- detecting overlap between matches and genes will work fine
                if start_coord < end_coord:
                    gene_list.append((start_coord, end_coord))
                else:
                    gene_list.append((end_coord, start_coord))
    
    return gene_list


def parse_args(sys_args: List[str]) -> argparse.Namespace:
    """
    Parses command line arguments.

    Args:
    - sys_args (List[str]): List of system arguments

    Returns:
    - parser.parse_args() (argparse.Namespace): Variables named after each argument, containing values parsed from
    command line
    """
    # TODO: Finish program, argument descriptions
    parser = argparse.ArgumentParser(sys_args, description="TODO")
    parser.add_argument("integration_tsv_dir", type=str, help="TODO")
    parser.add_argument("annotatation_tsv_dir", type=str, help="TODO")
    parser.add_argument("output_dir", type=str, help="TODO")
    parser.add_argument("--gene_count_threshold", type=int, help="TODO", default=DEFAULT_GENE_COUNT_THRESHOLD)
    parser.add_argument("--gene_coverage_threshold", type=float, help="TODO", default=DEFAULT_GENE_COVERAGE_THRESHOLD)
    parser.add_argument("--output_histogram", type=str, help="TODO", default="")
    parser.add_argument("--force", help="If output files already exist, overwrite them", action="store_true")

    return parser.parse_args()


def _main():
    # parse args from command line
    args = parse_args(sys.argv[1:])

    # path to directory containing output VIBES integration TSVs
    integration_tsv_dir = args.integration_tsv_dir

    # path to directory containing output VIBES viral gene annotation TSVs
    annotation_tsv_dir = args.annotation_tsv_dir

    # path to directory that will contain filtered output integration TSVs
    output_dir = args.output_dir

    # minimum number of genes required for a match to pass the filter
    minimum_gene_threshold = args.gene_count_threshold

    # how much of a gene must be contained in a match for it to be counted toward the filter threshold
    minimum_coverage = args.gene_coverage_threshold

    # path to output file to save histogram plot to
    hist_file_path = args.output_histogram

    # boolean that, when True, overwrites any output files should they exist
    force = args.force

    # dict where keys are the number of genes found in a match, and the number of matches with that many genes
    # are values
    histogram_dict = {}

    # get all file paths from integration directory
    integration_paths = list_files_in_dir(integration_tsv_dir)

    # get all file paths from gene annotation directory
    gene_annotation_paths = list_files_in_dir(annotation_tsv_dir)

    # get dict with gene position information for each query virus
    gene_dict = build_query_virus_gene_dict(gene_annotation_paths)

    for integration_path in integration_paths:
        with open(integration_path) as integration_file:
            # set output file path based on provided output directory and current integration file name
            output_file_path = f"{output_dir}/{Path(integration_path).stem}_filtered.tsv"
            output_file_contents = ""

            # some matches on different lines are part of the same integration, and share an ID. we want to count up
            # genes across all of these matches and include them as a group if they pass the threshold.
            prev_match_id = None
            prev_match_lines = ""
            prev_match_genes = 0

            for line in integration_file:
                # header lines start with #, so we grab it for the output file
                if line[0] == '#':
                    output_file_contents += line
                else:
                    match_tuple = generate_match_tuple(line)
                    genes_in_match = report_genes_in_match(match_tuple, gene_dict, minimum_coverage)

                    # check to see if the match id is the same as on the previous line. if so, continue to accumulate genes
                    # for this integration
                    if line[-1] == prev_match_id:
                        prev_match_genes += genes_in_match
                        prev_match_lines += line

                    # if the ids don't match, then we can see if the previous line(s) qualified:
                    else:
                        if prev_match_genes >= minimum_gene_threshold:
                            output_file_contents += prev_match_lines
                            # .setdefault() looks for the key in the dict, returning its value if the key exists.
                            # If not, it inserts the key with the second argument as a default value. So, here we
                            # add 1 to the value representing how many matches contained this many genes, adding
                            # 1 to 0 if it's the first case
                            histogram_dict[prev_match_genes] = histogram_dict.setdefault(prev_match_genes, 0) + 1

                    # either way, this line becomes the previous line
                    prev_match_id = line[-1]
                    prev_match_genes = genes_in_match
                    prev_match_lines = line

            # then, after exiting the loop, we know the previous line was the last list (and maybe more if it shared an
            # id with lines before it). We have to check if this last line passed the threshold:
            if prev_match_genes >= minimum_gene_threshold:
                output_file_contents += prev_match_lines
                # .setdefault() looks for the key in the dict, returning its value if the key exists.
                # If not, it inserts the key with the second argument as a default value. So, here we
                # add 1 to the value representing how many matches contained this many genes, adding
                # 1 to 0 if it's the first case
                histogram_dict[prev_match_genes] = histogram_dict.setdefault(prev_match_genes, 0) + 1

            # once we've iterated through all lines, we can write to the output file:
            write_lines_to_output(output_file_path, output_file_contents, force)

    plot_histogram(histogram_dict)