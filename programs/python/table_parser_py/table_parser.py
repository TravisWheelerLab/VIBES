import argparse
import os.path
import sys
import subprocess
from typing import TextIO
from typing import *
from os import path
from pathlib import Path
import json
import pandas as pd


INDENT_VAL = 4
TABLE_MODE = Literal["dfam", "tbl"]
JSON_MODE = Literal["integration", "annotation"]
STRAND = Literal["+", "-"]
# maximum allowable gap between two hits for them to be merge candidates
MAX_MERGE_DISTANCE = 5000



class QueryHit:
    def __init__(self, hit_name: str, acc_id: str, query_name: str, evalue: float, ali_st: int, ali_end: int, query_genome_path: str, hmm_st: int, hmm_end: int, hmm_len: int, strand: STRAND, verbose: bool, full_cutoff: int = 100, target_genome_len : int = None):
        # target = integration or (when annotating virus) protein, query = bacteria or (when annotating virus) virus
        self.hit_name = hit_name
        self.acc_id = acc_id
        self.query_name = query_name
        self.evalue = evalue
        self.ref_st = hmm_st
        self.ref_end = hmm_end
        self.query_genome_path = query_genome_path
        self.query_st = ali_st
        self.query_end = ali_end
        self.ref_len = hmm_len
        self.strand = strand
        self.verbose = verbose
        self.full_cutoff = full_cutoff

        if target_genome_len:
            self.target_genome_len = target_genome_len
        else:
            self.target_genome_len = self.get_genome_len()

        # TODO: Implement flanking att site attributes once detection implemented

    def get_seq_len(self) -> int:
        return abs(self.ref_st - self.ref_end) + 1

    def get_genome_len(self) -> int:
        seqstat_results = do_cmd(f"esl-seqstat {self.query_genome_path}", self.verbose)
        for line in seqstat_results.split("\n"):
            if "Total # residues" in line:
                line_list = line.split()
                return int(line_list[3])

    def get_percent_complete(self) -> float:
        return float(self.get_seq_len()) / self.ref_len

    def is_full_len(self) -> bool:
        return (self.ref_len - self.get_seq_len()) <= self.full_cutoff

    def detect_flanking_atts(self) -> Tuple[bool, List[str]]:
        # TODO: Implement this
        return None

    def to_tsv_line(self) -> str:
        return f"{self.hit_name}\t{self.acc_id}\t{self.query_name}\t{self.evalue}\t{self.is_full_len()}\t" \
               f"{self.query_st}\t{self.query_end}\t{self.ref_len}\t{path.basename(self.query_genome_path)}\t" \
               f"{self.ref_st}\t{self.ref_end}\t{self.target_genome_len}\t{self.strand}\n"


def strandedness(sense_pos: int, antisense_pos: int, strand: STRAND) -> int:
    if strand == "+":
        return sense_pos
    else:
        return antisense_pos


# This function tests whether current_hit is followed by next_hit on the reference viral genome
def order_preserved(current_hit: QueryHit, next_hit: QueryHit) -> bool:
    if current_hit.query_st < next_hit.query_st:
        return True
    else:
        return False


def should_merge(current_hit: QueryHit, next_hit: QueryHit) -> bool:
    current_ref_virus = current_hit.hit_name
    next_ref_virus = next_hit.hit_name

    if current_ref_virus == next_ref_virus:
        # get current end and next start positions, relative to bacterial genome. Depends on strand.
        current_bac_end = strandedness(current_hit.ref_end, current_hit.ref_st, current_hit.strand)
        next_bac_st = strandedness(next_hit.ref_st, next_hit.ref_end, next_hit.strand)

        bac_hit_distance = next_bac_st - current_bac_end

        if bac_hit_distance <= MAX_MERGE_DISTANCE and order_preserved(current_hit, next_hit):
            return True

    return False


# in some cases, one viral integration will be flagged by two or more hits. For instance, this may occur when regions of
# high homology are broken up by streches with little or no homology, encouraging the algorithm to call one integration
# as multiple unrelated hits. So we check for nearby, contiguous hits to the same reference viral genome. If spacing of
# the hits on the bacterial genome roughly matches the spacing of these same regions on the viral genome, we merge the
# hits.
def detect_merges(viral_hits: List[QueryHit]) -> List[QueryHit]:
    list_length = len(viral_hits)
    # for every hit in our list of detected viral sequences, compare to the next hit in the list. Merge the two if they
    # appear to be one viral integration incorrectly called as two separate hits
    for index, current_hit in enumerate(viral_hits):
        if (index + 1) < list_length:
            next_hit = viral_hits[index + 1]

            if should_merge(current_hit, next_hit):
                target_name = current_hit.hit_name
                query_name = current_hit.query_name
                # TODO: Eventually update to merge e-values
                evalue = max(current_hit.evalue, next_hit.evalue)
                # TODO: investigate prevalence of merge candidates with mismatched strandedness. However, we
                # TODO: don't anticipate this to be a major problem- nhmmer has no notion of the directionality of the
                # TODO: individual genes, so even antisense genes should be labeled as '+' if the instance in the
                # TODO: bacterial genome matches the direction of the gene in the reference viral genome. However, there
                # TODO: potentially may be cases where viral strains have 'flipped' genes, or something weird happened
                # TODO: to the integration post insertion. In the meantime, we just call strandedness to match hit 1
                strand = current_hit.strand
                bac_st = strandedness(current_hit.ref_st, next_hit.ref_st, strand)
                bac_end = strandedness(next_hit.ref_end, current_hit.ref_end, strand)
                bac_genome_path = current_hit.query_genome_path
                # if we've gotten this far, we already know that current_hit occurs before next_hit on the viral genome
                ref_vir_st = current_hit.query_st
                ref_vir_end = next_hit.query_end

                # TODO: create new merged hit, figure out how to get the list order sorted


def overwrite_check(file_path: str, force: bool) -> None:
    if path.isfile(file_path) and not force:
        raise FileExistsError(
            f"Output file {file_path} already exists- either move or delete this file or enable --force")


def write_occurrence_json(json_file: TextIO, query_hits: List[QueryHit]) -> None:
    # TODO: use Dict[str, List[int]] to store a list of ints for each virus present across the viral seqs. the list of ints
    # TODO: will be the occurrence count chart. This should easily plug into the json package
    occ_dict = {}

    for hit in query_hits:
        # each position in the list corresponds to a position in a reference target sequence (e.g. a viral
        # genome). the number at each position is how many times a nucleotide corresponding to that position has been
        # detected in a query sequence (e.g. a bacterial genome)
        if hit.hit_name not in occ_dict.keys():
            occ_dict[hit.hit_name] = [0] * hit.ref_len

        ref_start_index = hit.ref_st - 1 # offset by 1 because genomes are 1-indexed, lists are 0-indexed
        ref_end_index = hit.ref_end # don't offset by 1 because range() excludes end

        for index in range(ref_start_index, ref_end_index):
            occ_dict[hit.hit_name][index] += 1

    json_file.write(json.dumps(occ_dict, indent=INDENT_VAL))


def write_occurrence_json_from_path(json_path: str, query_hits: List[QueryHit], force: bool) -> None:
    overwrite_check(json_path, force)

    with open(json_path, "w") as json_file:
        write_occurrence_json(json_file, query_hits)
        json_file.close()


def get_genome_len(genome_path: str, verbose: bool) -> int:
    seqstat_results = do_cmd(f"esl-seqstat {genome_path}", verbose)
    for line in seqstat_results.split("\n"):
        if "Total # residues" in line:
            line_list = line.split()
            return int(line_list[3])


def write_annotation_json(json_file: TextIO, query_hits: List[QueryHit], protein_annotations: Optional[Dict[str, str]], occurrence_dict: Dict[str, List[int]], genome_path: str, verbose: bool) -> None:
    phage_name = Path(genome_path).with_suffix('').stem

    if phage_name in occurrence_dict.keys():
        annotation_json = {"prophageName": phage_name,
                           "occurrences": occurrence_dict[phage_name]}
    else:
        # if this function in being run, we expect the genome to be the reference phage genome
        genome_len = get_genome_len(genome_path, verbose)
        annotation_json = {"prophageName": phage_name,
                           "occurrences": [0] * genome_len}

    prot_anno_list = []

    for hit in query_hits:
        hit_dict = {"name": hit.hit_name,
                    "ID": hit.acc_id,
                    "start": hit.query_st,
                    "end": hit.query_end,
                    "evalue": hit.evalue,
                    "desc": protein_annotations[hit.hit_name]
                    }
        prot_anno_list.append(hit_dict)

    annotation_json["annotations"] = prot_anno_list

    json_file.write(json.dumps(annotation_json))

def write_annotation_json_from_path(json_path: str, query_hits: List[QueryHit], protein_annotations: Optional[Dict[str, str]], occurrence_dict: Dict[str, List[int]], genome_path: str, force: bool, verbose: bool) -> None:
    overwrite_check(json_path, force)

    with open(json_path, "w") as json_file:
        write_annotation_json(json_file, query_hits, protein_annotations, occurrence_dict, genome_path, verbose)
        json_file.close()


def write_tsv(tsv: TextIO, seq_list: List[QueryHit]) -> None:
    # write header line first
    headers = ["Hit Name", "Accession", "Query Sequence Name", "E-Value", "Full Length",
               "Ref Seq Start", "Ref Seq End", "Ref Seq Length",
               "Bacterial Genome Name", "Match Start on Query Seq",
               "Match End on Query Seq", "Query Genome Length", "Strand\n"]
    tsv.write("\t".join(headers))
    for viral_seq in seq_list:
        tsv.write(viral_seq.to_tsv_line())


def write_tsv_from_path(tsv_path: str, seq_list: List[QueryHit], force: bool) -> None:
    overwrite_check(tsv_path, force)

    with open(tsv_path, "w") as tsv:
        write_tsv(tsv, seq_list)


def parse_dfam(dfam_file: TextIO, genome_path: str, max_eval: float, verbose) -> List[QueryHit]:
    hit_list = []
    for line_num, line in enumerate(dfam_file, 0):
        if line[0] == "#":
            pass
        else:
            line_list = line.split()

            target_name = line_list[0]
            acc_id = str(line_list[1])
            query_name = line_list[2]
            evalue = float(line_list[4])
            hmm_st = int(line_list[6])
            hmm_en = int(line_list[7])
            strand = line_list[8]
            ali_st = int(line_list[9])
            ali_en = int(line_list[10])
            hmm_len = int(line_list[13])

            if evalue <= max_eval:
                hit_list.append(QueryHit(target_name, acc_id, query_name, evalue, ali_st, ali_en, genome_path, hmm_st, hmm_en, hmm_len, strand, verbose))
            else:
                if verbose:
                    print(f"Excluding line {line_num}: e-value of {evalue} failed to pass maximum e-value threshold of {max_eval}")

    return hit_list

def parse_tbl(tbl_file: TextIO, genome_path: str, max_eval: float, verbose: bool) -> List[QueryHit]:
    hit_list = []
    for line_num, line in enumerate(tbl_file, 0):
        if line[0] == "#":
            pass
        else:
            line_list = line.split()

            # FraHMMER and nhmmscan disagree over what is the target and what is the query- we use nhmmscan's notation
            target_name = str(line_list[0])
            acc_id = str(line_list[1])
            query_name = str(line_list[2])
            hmm_len = int(line_list[4])
            hmm_st = int(line_list[5])
            hmm_en = int(line_list[6])
            ali_st = int(line_list[8])
            ali_en = int(line_list[9])
            evalue = float(line_list[12])

            strand = "+"
            if hmm_en < hmm_st:
                strand = "-"

            if evalue <= max_eval:
                hit_list.append(QueryHit(query_name, acc_id, target_name, evalue, ali_st, ali_en, genome_path, hmm_st, hmm_en, hmm_len, strand, verbose))
            else:
                if verbose:
                    print(f"Excluding line {line_num}: e-value of {evalue} failed to pass maximum e-value threshold of {max_eval}")

    return hit_list


def parse_table_from_path(table_path: str, genome_path: str, max_eval: float, table_mode: TABLE_MODE,  verbose: bool) -> List[QueryHit]:
    with open(table_path) as table_file:
        if verbose:
            print(f"Opening {table_path}...")

        if table_mode == "dfam":
            return parse_dfam(table_file, genome_path, max_eval, verbose)
        elif table_mode == "tbl":
            return parse_tbl(table_file, genome_path, max_eval, verbose)
        else:
            raise ValueError("table_type must be either dfam or tbl")


def parse_protein_annotation_from_path(anno_tsv_path: str, verbose: bool) -> Dict[str, str]:
    if verbose:
        print(f"Opening {anno_tsv_path} with Pandas...")

    anno_dict = {}

    with open(anno_tsv_path, "r") as anno_file:
        for line in anno_file:
            phrog, desc = line.split("\t")
            anno_dict[phrog.rstrip()] = desc.rstrip()

    return anno_dict


def load_json(json_path: str, verbose: bool) -> Dict[str, List[int]]:
    if verbose:
        print(f"Opening {json_path}...")

    with open(json_path, "r") as json_file:
        return json.load(json_file)


def parse_args(sys_args: list) -> argparse.Namespace:
    default_eval = 1e-5
    parser = argparse.ArgumentParser(sys_args, description="Parses input .dfam or .tbl (from FraHMMER --tblout) file, extracting information about query hits detected in target genome")
    parser.add_argument("table_path", type=str, help="Path to input .dfam or .tbl file made up of query hits")
    parser.add_argument("genome_path", type=str, help="Path to target genome in .fasta format")
    parser.add_argument("output_tsv_path", type=str, help="Path to output .tsv file")
    parser.add_argument("output_json_path", type=str, help="Path to output .json file containing information on "
                                                           "nucleotide occurrence counts")
    parser.add_argument("table_type", type=str, choices=["dfam","tbl"], default="dfam", help="Which type of table is being supplied as input, which must be dfam or tbl (default dfam)")
    parser.add_argument("json_type", type=str, choices=["integration","annotation"], default="integration", help="Sets which type of .json will be produced: integration, which counts how many times a position in a user-supplied viral genome has been detected in a specific bacterial genome, or annotation, which produces .json describing a reference viral genome annotated with viral proteins. Default is integration")
    parser.add_argument("--occurrence_json", type=str, default="", help="Path to summed occurrence .json. This file contains a field for each virus found integrated at least once across all bacteria, which contains a per-position count of how many times each position has been detected in an integration. Necessary to produce annotation .jsons")
    parser.add_argument("--annotation_tsv", type=str, default="", help="Path to .tsv file containing information about the function of viral proteins used to annotate user-supplied viruses")
    parser.add_argument("--max_evalue", type=float, default=default_eval, help=f"Maximum allowed sequence e-value (must be >= 0). Default is {default_eval}")
    parser.add_argument("--verbose", action="store_true", help="Print additional information useful for debugging")
    parser.add_argument("--force", action="store_true", help="If output file already exists, overwrite it")

    return parser.parse_args()


def do_cmd(cmd: str, verbose: bool) -> str:
    if verbose:
        print(f"Running command: {cmd}")

    return subprocess.getoutput(cmd)


def _main():
    args = parse_args(sys.argv[1:])
    table_path = args.table_path
    genome_path = args.genome_path
    tsv_path = args.output_tsv_path
    json_path = args.output_json_path
    table_mode = args.table_type
    json_mode = args.json_type
    occ_json_path = args.occurrence_json
    protein_annotation_path = args.annotation_tsv
    max_eval = args.max_evalue
    verbose = args.verbose
    force = args.force

    # check that inputs are legal
    if max_eval < 0:
        raise ValueError("--max_evalue must be used with an argument greater than or equal to 0")

    # integration mode
    if json_mode == "integration":
        query_hits = parse_table_from_path(table_path, genome_path, max_eval, table_mode, verbose)
        write_tsv_from_path(tsv_path, query_hits, force)
        write_occurrence_json_from_path(json_path, query_hits, force)

    # annotation mode
    elif json_mode == "annotation":
        protein_annotations = None

        if protein_annotation_path:
            protein_annotations = parse_protein_annotation_from_path(protein_annotation_path, verbose)

        query_hits = parse_table_from_path(table_path, genome_path, max_eval, table_mode, verbose)
        occurrence_dict = load_json(occ_json_path, verbose)
        write_tsv_from_path(tsv_path, query_hits, force)
        write_annotation_json_from_path(json_path, query_hits, protein_annotations, occurrence_dict, genome_path, force, verbose)






if __name__ == "__main__":
    _main()
