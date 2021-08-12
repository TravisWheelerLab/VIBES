import argparse
import sys
import subprocess
from typing import *
from os import path
from os import remove


STRAND = Literal["+", "-"]


class ViralSeq:
    def __init__(self, name: str, evalue: float, bac_st: int, bac_end: int, bac_genome_path: str, ref_vir_st: int, ref_vir_end: int, ref_vir_len: int, strand: STRAND, full_cutoff: int = 100):
        self.name = name
        self.evalue = evalue
        self.bac_st = bac_st
        self.bac_end = bac_end
        self.bac_genome_path = bac_genome_path
        self.ref_vir_st = ref_vir_st
        self.ref_vir_end = ref_vir_end
        self.ref_vir_len = ref_vir_len
        self.strand = strand
        self.full_cutoff = full_cutoff
        self.len = self.__buildLen__()
        self.bac_genome_len = self.__fetch_bac_len__()
        self.flanking_att_site = self.__detect_flanking_atts__()

    def get_vir_seq_len(self) -> int:
        return abs(self.gn_st - self.gn_end) + 1

    def get_bac_len(self) -> int:
        seqstat_results = do_cmd(f"esl-seqstat {self.bac_genome_path}").stdout
        for line in seqstat_results:
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

def parse_table(dfam_file: TextIO, genome_path: str) -> List[ViralSeq]:
    seq_list = []
    for line in dfam_file:
        line_list = line.split()
        name = line_list[0]
        evalue = line_list[4]
        hmm_st = line_list[6]
        hmm_en = line_list[7]
        strand = line_list[8]
        ali_st = line_list[9]
        ali_en = line_list[10]
        ref_vir_len = line_list[13]
        # target name 0, e-value 4, start on target phage (hmmst) 6, end on target phage (hmmen) 7, strand 8, query bacteria start (alist) 9, query bacteria end (alien) 10, phage len (modlen) 13

        seq_list.append(ViralSeq(name, evalue, ali_st, ali_en, genome_path, hmm_st, hmm_en, ref_vir_len, strand))

    return seq_list


def parse_table_from_path(dfam_path: str, genome_path: str) -> List[ViralSeq]:
    with open(dfam_path) as dfam_file:
        parse_table(dfam_file)


def parse_args(sys_args: list) -> argparse.Namespace:
    default_eval = 1e-5
    parser = argparse.ArgumentParser(sys_args, description="Parses input dfam file, extracting information about viral sequence in the bacterial genome")
    parser.add_argument("dfam_path", type=str, help="Path to input .dfam file")
    parser.add_argument("genome_path", type=str, help="Path to input genome in .fasta format")
    parser.add_argument("output_tsv_path", type=str, help="Path to output .tsv file")
    parser.add_argument("output_json_path", type=str, help="Path to output .json file containing information on nucleotide occurrence counts")
    parser.add_argument("--max_evalue", type=float, default=default_eval, help=f"Maximum allowable sequence e-value. Default is {default_eval}")
    parser.add_argument("--verbose", action="store_true", help="Print additional information useful for debugging")
    parser.add_argument("--force", action="store_true", help="If output file already exists, overwrite it")

    return parser.parse_args()


def do_cmd(cmd: List[str], verbose: bool) -> subprocess.CompletedProcess:
    # double check that all elements are strings
    cmd = [str(e) for e in cmd]
    if verbose:
        verbose_cmd = " ".join(cmd)
        print(f"Running command: {verbose_cmd}")

    return subprocess.run(cmd)


def _main():
    args = parse_args(sys.argv[1:])
    hmm_path = args.hmm_path
    genome_path = args.genome_path
    dfam_dir = args.output_table_dir
    cpu_count = args.cpu
    verbose = args.verbose
    force = args.force

    # check that inputs are legal
    if cpu_count < 0:
        raise ValueError("--cpu must be used with an argument greater than or equal to 0")

    if not path.isfile(hmm_path):
        raise FileNotFoundError(f"No such file: {hmm_path}")

    if not path.isfile(genome_path):
        raise FileNotFoundError(f"No such file: {genome_path}")

    if not path.isdir(dfam_dir):
        raise NotADirectoryError(f"No such folder: {dfam_dir}")


if __name__ == "__main__":
    _main()