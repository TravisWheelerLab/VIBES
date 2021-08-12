import argparse
import sys
import subprocess
from typing import TextIO
from typing import *
from os import path
from os import remove


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
        self.flanking_att_site = self.detect_flanking_atts()

    def get_vir_seq_len(self) -> int:
        return abs(self.gn_st - self.gn_end) + 1

    def get_bac_len(self) -> int:
        seqstat_results = do_cmd(f"esl-seqstat {self.bac_genome_path}", self.verbose)
        for line in seqstat_results.split("\n"):
            if "Total # residues" in line:
                line_list = line.split()
                return line_list[3]

    def get_percent_complete(self) -> float:
        return float(self.len) / self.ref_vir_len

    def is_full_len(self) -> bool:
        return (self.ref_vir_len - self.len) <= self.full_cutoff

    def detect_flanking_atts(self):
        # TODO: use simple smith-waterman alignment to search for off-diagonal hits at the ends of viral genome
        pass

def parse_table(dfam_file: TextIO, genome_path: str, verbose) -> List[ViralSeq]:
    seq_list = []
    for line in dfam_file:
        if line[0] == "#":
            pass
        else:
            line_list = line.split()
            name = line_list[0]
            evalue = line_list[4]
            hmm_st = line_list[6]
            hmm_en = line_list[7]
            strand = line_list[8]
            ali_st = line_list[9]
            ali_en = line_list[10]
            ref_vir_len = line_list[13]

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
    parser.add_argument("output_json_path", type=str, help="Path to output .json file containing information on nucleotide occurrence counts")
    parser.add_argument("--max_evalue", type=float, default=default_eval, help=f"Maximum allowable sequence e-value (must be >= 0). Default is {default_eval}")
    parser.add_argument("--verbose", action="store_true", help="Print additional information useful for debugging")
    parser.add_argument("--force", action="store_true", help="If output file already exists, overwrite it")

    return parser.parse_args()


def do_cmd(cmd: str, verbose: bool) -> subprocess.CompletedProcess:
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
    print(viral_seqs[0].bac_genome_len)


if __name__ == "__main__":
    _main()