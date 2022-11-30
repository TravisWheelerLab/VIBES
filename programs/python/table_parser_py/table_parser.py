import argparse
import sys
import subprocess
from typing import TextIO
from typing import *
from os import path
from os import remove
import json


STRAND = Literal["+", "-"]
# maximum allowable gap between two hits for them to be merge candidates
MAX_MERGE_DISTANCE = 5000



class ViralSeq:
    def __init__(self, name: str, query_name: str, evalue: float, bac_st: int, bac_end: int, bac_genome_path: str, ref_vir_st: int, ref_vir_end: int, ref_vir_len: int, strand: STRAND, full_cutoff: int = 100, verbose: bool = False):
        self.name = name
        self.query_name = query_name
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
        return f"{self.name}\t{self.query_name}\t{self.evalue}\t{self.is_full_len()}\t{self.ref_vir_st}\t" \
               f"{self.ref_vir_end}\t{self.ref_vir_len}\t{self.is_flanked}\t{path.basename(self.bac_genome_path)}\t" \
               f"{self.bac_st}\t{self.bac_end}\t{self.bac_genome_len}\t{self.strand}\n"

def strandedness(sense_pos: int, antisense_pos: int, strand: STRAND) -> int:
    if strand == "+":
        return sense_pos
    else:
        return antisense_pos

# This function tests whether current_hit is followed by next_hit on the reference viral genome
def vir_order_preserved(current_hit: ViralSeq, next_hit: ViralSeq) -> bool:
    if current_hit.ref_vir_st < next_hit.ref_vir_st:
        return True
    else:
        return False

def should_merge(current_hit: ViralSeq, next_hit: ViralSeq) -> bool:
    current_ref_virus = current_hit.name
    next_ref_virus = next_hit.name

    if current_ref_virus == next_ref_virus:
        # get current end and next start positions, relative to bacterial genome. Depends on strand.
        current_bac_end = strandedness(current_hit.bac_end, current_hit.bac_st, current_hit.strand)
        next_bac_st = strandedness(next_hit.bac_st, next_hit.bac_end, next_hit.strand)

        bac_hit_distance = next_bac_st - current_bac_end

        if bac_hit_distance <= MAX_MERGE_DISTANCE and vir_order_preserved(current_hit, next_hit):
            return True

    return False


# in some cases, one viral integration will be flagged by two or more hits. For instance, this may occur when regions of
# high homology are broken up by streches with little or no homology, encouraging the algorithm to call one integration
# as multiple unrelated hits. So we check for nearby, contiguous hits to the same reference viral genome. If spacing of
# the hits on the bacterial genome roughly matches the spacing of these same regions on the viral genome, we merge the
# hits.
def detect_merges(viral_hits: List[ViralSeq]) -> List[ViralSeq]:
    list_length = len(viral_hits)
    # for every hit in our list of detected viral sequences, compare to the next hit in the list. Merge the two if they
    # appear to be one viral integration incorrectly called as two separate hits
    for index, current_hit in enumerate(viral_hits):
        if (index + 1) < list_length:
            next_hit = viral_hits[index + 1]

            if should_merge(current_hit, next_hit):
                target_name = current_hit.name
                query_name = current_hit.query_name
                # TODO: Discuss how to merge evalues with Travis
                #evalue = TODO
                #strand = TODO: investigate prevalence of merge candidates with mismatched strandedness. However, we
                # TODO: don't anticipate this to be a major problem- nhmmer has no notion of the directionality of the
                # TODO: individual genes, so even antisense genes should be labeled as '+' if the instance in the
                # TODO: bacterial genome matches the direction of the gene in the reference viral genome. However, there
                # TODO: potentially may be cases where viral strains have 'flipped' genes, or something weird happened
                # TODO: to the integration post insertion. In the meantime, we just call strandedness to match hit 1
                strand = current_hit.strand
                bac_st = strandedness(current_hit.bac_st, next_hit.bac_st, strand)
                bac_end = strandedness(next_hit.bac_end, current_hit.bac_end, strand)
                bac_genome_path = current_hit.bac_genome_path
                # if we've gotten this far, we already know that current_hit occurs before next_hit on the viral genome
                ref_vir_st = current_hit.ref_vir_st
                ref_vir_end = next_hit.ref_vir_end



def populate_occurrence_json(json_file: TextIO, viral_seqs: List[ViralSeq]) -> None:
    # TODO: use Dict[str, List[int]] to store a list of ints for each virus present across the viral seqs. the list of ints
    # TODO: will be the occurrence count chart. This should easily plug into the json package
    occ_dict = {}
    INDENT_VAL = 4

    for seq in viral_seqs:
        # each position in the list corresponds to a position in the reference viral genome. the number at each
        # position is how many times we've seen that position as part of an identified viral sequence
        if seq.name not in occ_dict.keys():
            occ_dict[seq.name] = [0] * seq.ref_vir_len

        ref_start_index = seq.ref_vir_st - 1 # offset by 1 because genomes are 1-indexed, lists are 0-indexed
        ref_end_index = seq.ref_vir_end # don't offset by 1 because range() excludes end

        for index in range(ref_start_index, ref_end_index):
            occ_dict[seq.name][index] += 1

    json_file.write(json.dumps(occ_dict, indent=INDENT_VAL))


def populate_occurrence_json_from_path(json_path: str, viral_seqs: List[ViralSeq], force: bool) -> None:
    if path.isfile(json_path) and not force:
        raise FileExistsError(
            f"Output file {json_path} already exists- either move or delete this file or enable --force")
    else:
        with open(json_path, "w") as json_file:
            populate_occurrence_json(json_file, viral_seqs)
            json_file.close()


def populate_tsv(tsv: TextIO, seq_list: List[ViralSeq]) -> None:
    # write header line first
    headers = ["Name", "Query Sequence Name", "E-Value", "Is Full Length Insertion",
               "Match Start Position on Viral Genome", "Match End Position on Viral Genome", "Viral Genome Length",
               "Flanking Att Sites", "Bacterial Genome Name", "Match Start on Query Seq",
               "Match End on Query Seq", "Query Genome Length", "Strand\n"]
    tsv.write("\t".join(headers))
    for viral_seq in seq_list:
        tsv.write(viral_seq.to_tsv_line())


def populate_tsv_from_path(tsv_path: str, seq_list: List[ViralSeq], force: bool) -> None:
    if path.isfile(tsv_path) and not force:
        raise FileExistsError(
            f"Output file {tsv_path} already exists- either move or delete this file or enable --force")
    else:
        with open(tsv_path, "w") as tsv:
            populate_tsv(tsv, seq_list)


def parse_table(dfam_file: TextIO, genome_path: str, max_eval: float, verbose) -> List[ViralSeq]:
    seq_list = []
    for line_num, line in enumerate(dfam_file, 0):
        if line[0] == "#":
            pass
        else:
            line_list = line.split()
            target_name = line_list[0]
            query_name = line_list[2]
            evalue = float(line_list[4])
            hmm_st = int(line_list[6])
            hmm_en = int(line_list[7])
            strand = line_list[8]
            ali_st = int(line_list[9])
            ali_en = int(line_list[10])
            ref_vir_len = int(line_list[13])

            if evalue <= max_eval:
                seq_list.append(ViralSeq(target_name, query_name, evalue, ali_st, ali_en, genome_path, hmm_st, hmm_en, ref_vir_len, strand, verbose))
            else:
                if verbose:
                    print(f"Excluding line {line_num}: e-value of {evalue} failed to pass maximum e-value threshold of {max_eval}")

    return seq_list


def parse_table_from_path(dfam_path: str, genome_path: str, max_eval: float, verbose) -> List[ViralSeq]:
    with open(dfam_path) as dfam_file:
        if verbose:
            print(f"Opening {dfam_path}...")
        return parse_table(dfam_file, genome_path, max_eval, verbose)


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

    viral_seqs = parse_table_from_path(dfam_path, genome_path, max_eval, verbose)
    populate_tsv_from_path(tsv_path, viral_seqs, force)
    populate_occurrence_json_from_path(json_path, viral_seqs, force)


if __name__ == "__main__":
    _main()
