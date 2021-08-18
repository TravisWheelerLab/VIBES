import argparse
import sys
import subprocess
from typing import TextIO
from typing import *
from os import path
from os import remove
import json


STRAND = Literal["+", "-"]


class ViralSeq:
    def __init__(self, name: str, evalue: float, bac_st: int, bac_end: int, bac_genome_path: str, ref_vir_st: int, ref_vir_end: int, ref_vir_len: int, strand: STRAND, full_cutoff: int = 100, verbose: bool = False):
        self.name = name
        self.evalue = evalue
        self.bac_st = bac_st
        self.bac_end = bac_end
        self.bac_genome_path = bac_genome_path
        self.ref_vir_st = ref_vir_st
        self.ref_vir_end = ref_vir_end
        self.ref_vir_len = ref_vir_len
        self.strand = strand
        self.verbose = verbose
        self.full_cutoff = full_cutoff
        self.bac_genome_len = self.get_bac_len()
        self. is_flanked, self.flanking_att_site = self.detect_flanking_atts()

    def get_seq_len(self) -> int:
        return abs(self.bac_st - self.bac_end) + 1

    def get_bac_len(self) -> int:
        seqstat_results = do_cmd(f"esl-seqstat {self.bac_genome_path}", self.verbose)
        for line in seqstat_results.split("\n"):
            if "Total # residues" in line:
                line_list = line.split()
                return int(line_list[3])

    def get_percent_complete(self) -> float:
        return float(self.get_seq_len()) / self.ref_vir_len

    def is_full_len(self) -> bool:
        return (self.ref_vir_len - self.get_seq_len()) <= self.full_cutoff

    def detect_flanking_atts(self) -> Tuple[bool, List[str]]:
        # TODO: use simple smith-waterman alignment to search for off-diagonal hits at the ends of viral genome
        # in the meantime, simply return as negative
        return False, ["",""]

    def to_tsv_line(self) -> str:
        return f"{self.name}\t{self.evalue}\t{self.is_full_len()}\t{self.ref_vir_st}\t{self.ref_vir_end}\t" \
               f"{self.ref_vir_len}\t{self.is_flanked}\t{path.basename(self.bac_genome_path)}\t{self.bac_st}\t" \
               f"{self.bac_end}\t{self.bac_genome_len}\t{self.strand}\n"


def populate_occurrence_json(json_file: TextIO, viral_seqs: List[ViralSeq]):
    # TODO: use Dict[str, List[int]] to store a list of ints for each virus present across the viral seqs. the list of ints
    # will be the occurrence count chart. This should easily plug into the json package
    occ_dict = {}
    INDENT_VAL = 4

    for seq in viral_seqs:
        # each position in the list corresponds to a position in the reference viral genome. the number at each
        # position is how many times we've seen that position as part of an identified viral sequence
        if seq.name not in occ_dict:
            occ_dict[seq.name] = [0] * seq.ref_vir_len

        ref_start_index = seq.ref_vir_st - 1 # offset by 1 because genomes are 1-indexed, lists are 0-indexed
        ref_end_index = seq.ref_vir_end # don't offset by 1 because range() excludes end

        for index in range(ref_start_index, ref_end_index):
            occ_dict[seq.name][index] += 1

        json_file.write(json.dumps(occ_dict, indent=INDENT_VAL))


def populate_occurrence_json_from_path(json_path: str, viral_seqs: List[ViralSeq], force: bool):
    if path.isfile(json_path) and not force:
        raise FileExistsError(
            f"Output file {json_path} already exists- either move or delete this file or enable --force")
    else:
        with open(json_path, "w") as json_file:
            populate_occurrence_json(json_file, viral_seqs)
            json_file.close()


def populate_tsv(tsv: TextIO, seq_list: List[ViralSeq]):
    # write header line first
    headers = ["Name", "E-Value", "Is Full Length Insertion", "Match Start Position on Viral Genome",
               "Match End Position on Viral Genome", "Viral Genome Length", "Flanking Att Sites",
               "Bacterial Genome Name", "Match Start Position on Bacterial Genome",
               "Match End Position on Bacterial Genome", "Bacterial Genome Length", "Strand\n"]
    tsv.write("\t".join(headers))
    for viral_seq in seq_list:
        tsv.write(viral_seq.to_tsv_line())


def populate_tsv_from_path(tsv_path: str, seq_list: List[ViralSeq], force: bool):
    if path.isfile(tsv_path) and not force:
        raise FileExistsError(
            f"Output file {tsv_path} already exists- either move or delete this file or enable --force")
    else:
        with open(tsv_path, "w") as tsv:
            populate_tsv(tsv, seq_list)


def parse_table(dfam_file: TextIO, genome_path: str, verbose) -> List[ViralSeq]:
    seq_list = []
    for line in dfam_file:
        if line[0] == "#":
            pass
        else:
            line_list = line.split()
            name = line_list[0]
            evalue = float(line_list[4])
            hmm_st = int(line_list[6])
            hmm_en = int(line_list[7])
            strand = line_list[8]
            ali_st = int(line_list[9])
            ali_en = int(line_list[10])
            ref_vir_len = int(line_list[13])

            seq_list.append(ViralSeq(name, evalue, ali_st, ali_en, genome_path, hmm_st, hmm_en, ref_vir_len, strand, verbose))

    return seq_list


def parse_table_from_path(dfam_path: str, genome_path: str, verbose) -> List[ViralSeq]:
    with open(dfam_path) as dfam_file:
        return parse_table(dfam_file, genome_path, verbose)


def parse_args(sys_args: list) -> argparse.Namespace:
    default_eval = 1e-5
    parser = argparse.ArgumentParser(sys_args, description="Parses input dfam file, extracting information about viral sequence in the bacterial genome")
    parser.add_argument("dfam_path", type=str, help="Path to input .dfam file")
    parser.add_argument("genome_path", type=str, help="Path to input genome in .fasta format")
    parser.add_argument("output_tsv_path", type=str, help="Path to output .tsv file")
    parser.add_argument("output_json_path", type=str, help="Path to output .json file containing information on "
                                                           "nucleotide occurrence counts")
    parser.add_argument("--max_evalue", type=float, default=default_eval, help=f"Maximum allowable sequence e-value (must be >= 0). Default is {default_eval}")
    parser.add_argument("--verbose", action="store_true", help="Print additional information useful for debugging")
    parser.add_argument("--force", action="store_true", help="If output file already exists, overwrite it")

    return parser.parse_args()


def do_cmd(cmd: str, verbose: bool) -> str:
    if verbose:
        print(f"Running command: {cmd}")

    return subprocess.getoutput(cmd)


def _main():
    args = parse_args(sys.argv[1:])
    dfam_path = args.dfam_path
    genome_path = args.genome_path
    tsv_path = args.output_tsv_path
    json_path = args.output_json_path
    max_eval = args.max_evalue
    verbose = args.verbose
    force = args.force

    # check that inputs are legal
    if max_eval < 0:
        raise ValueError("--max_evalue must be used with an argument greater than or equal to 0")

    viral_seqs = parse_table_from_path(dfam_path, genome_path, verbose)
    populate_tsv_from_path(tsv_path, viral_seqs, force)
    populate_occurrence_json_from_path(json_path, viral_seqs, force)


if __name__ == "__main__":
    _main()