import argparse
import sys
import subprocess
from typing import TextIO
from typing import *
from os import path
from pathlib import Path
import json
# Dependency: esl-seqstat in some circumstances

INDENT_VAL = 4
TABLE_MODE = Literal["dfam", "tbl"]
JSON_MODE = Literal["integration", "annotation"]
STRAND = Literal["+", "-"]

# maximum percentage of reference genome length we allow to occur between merge candidates on reference (viral or
# protein) seq
REFERENCE_MAX_GAP_FACTOR = .25
# TODO: Make this a command line option


# TODO: Document this class and its quirks
class QueryHit:
    def __init__(self, hit_name: str, acc_id: str, query_name: str, evalue: float, ali_st: int, ali_end: int,
                 query_genome_name: str, hmm_st: int, hmm_end: int, hmm_len: int, strand: STRAND, verbose: bool,
                 query_genome_len : int = None, description: str = "-"):
        # target = integration or (when annotating virus) protein, query = bacteria or (when annotating virus) virus
        self.hit_name = hit_name
        self.acc_id = acc_id
        self.query_name = query_name
        self.evalue = evalue
        self.ref_st = hmm_st
        self.ref_end = hmm_end
        self.query_genome_name = query_genome_name
        self.query_st = ali_st
        self.query_end = ali_end
        self.ref_len = hmm_len
        self.strand = strand
        self.verbose = verbose
        self.full_length = None
        self.description = description
        self.integration_id = "" # Document quirk here

        if query_genome_len:
            self.query_genome_len = query_genome_len
        else:
            self.query_genome_len = self.get_genome_len()

    def get_seq_len_on_ref(self) -> int:
        # these positions are 1-indexed, so we have to add one
        # (if you start at position 1 and go to position 5, the sequence is 5 positions long, not 4)
        return abs(self.ref_end - self.ref_st) + 1

    def get_genome_len(self) -> int:
        seqstat_results = do_cmd(f"esl-seqstat {self.query_genome_name}", self.verbose)
        for line in seqstat_results.split("\n"):
            if "Total # residues" in line:
                line_list = line.split()
                return int(line_list[3])

    def get_percent_complete(self) -> float:
        return float(self.get_seq_len_on_ref()) / self.ref_len

# TODO: Look up deprecation annotation in Python

    def to_tsv_line(self) -> str:
        return f"{self.hit_name}\t{self.description}\t{self.acc_id}\t{self.evalue}\t{self.full_length}\t" \
               f"{self.ref_st}\t{self.ref_end}\t{self.ref_len}\t{path.basename(self.query_genome_name)}\t" \
               f"{self.query_name}\t{self.query_st}\t{self.query_end}\t{self.query_genome_len}\t{self.strand}\t" \
               f"{self.integration_id}\n"


def find_position_for_strand_type(sense_position: int, antisense_position: int, strand: STRAND) -> int:
    if strand == "+":
        return sense_position
    else:
        return antisense_position


# This function tests whether current_hit is followed by next_hit on the reference viral genome
def ref_order_preserved(current_hit: QueryHit, prev_hit: QueryHit, overlap_tolerance_percent) -> bool:
    # if strand is +, then prev_hit occurs before current_hit on both the reference viral and bacterial genomes.
    # if strand is -, then current_hit occurs before prev_hit on the viral genome, but prev_hit is still first on the
    # bacterial genome (which is why it's the previous hit)

    # generally, we expect that the end of the preceeding hit on the viral genome comes before the start of the next hit
    # in which case overlap will be negative
    # if a significant mutation has occurred, then the second hit's viral genome start position might come slightly
    # before the first hit's end position, so we allow for a bit of overlap
    overlap = find_position_for_strand_type(prev_hit.ref_end, current_hit.ref_end, prev_hit.strand) \
            - find_position_for_strand_type(current_hit.ref_st, prev_hit.ref_st, prev_hit.strand)
    max_overlap = current_hit.ref_len * overlap_tolerance_percent

    return overlap <= max_overlap


def names_and_strand_match(prev_hit: QueryHit, current_hit: QueryHit) -> bool:
    return current_hit.hit_name == prev_hit.hit_name and current_hit.strand == prev_hit.strand and \
        current_hit.query_name == prev_hit.query_name


def same_integration(prev_hit: QueryHit, current_hit: QueryHit, overlap_tolerance_percent: float) -> bool:
    strand = current_hit.strand

    # if both match to same reference virus
    # TODO: make this evaluation a function call
    if names_and_strand_match(prev_hit, current_hit):
        # get where previous hit ends and the current hit starts, from the perspective of the bacterial genome
        prev_bac_end = find_position_for_strand_type(prev_hit.query_end, prev_hit.query_st, strand)
        current_bac_st = find_position_for_strand_type(current_hit.query_st, current_hit.query_end, strand)
        # if strand is +, then prev_hit should have the lower reference start position. if strand is -, then current hit
        # will have the lower reference start position (since it's reversed from the perspective of the bacterial genome)
        # note: this will only be the case if fragments are actually part of a larger insertion. May not hold if they're
        # unrelated
        lower_ref_end = find_position_for_strand_type(prev_hit.ref_end, current_hit.ref_end, strand)
        higher_ref_st = find_position_for_strand_type(current_hit.ref_st, prev_hit.ref_st, strand)

        max_hit_gap = current_hit.ref_len * REFERENCE_MAX_GAP_FACTOR
        bac_hit_gap = current_bac_st - prev_bac_end
        # if the gaps are small enough* and hits are in the right order:
        # *some unrelated hits won't obey the logic above, potentially resulting in negative gap values
        # TODO: ake this a function too
        if 0 <= bac_hit_gap <= max_hit_gap and ref_order_preserved(current_hit, prev_hit, overlap_tolerance_percent):
            return True

    return False


# in some cases, one viral integration will be flagged by two or more hits. For instance, this may occur when regions of
# high homology are broken up by streches with little or no homology, encouraging the algorithm to call one integration
# as multiple unrelated hits. So we check for nearby, contiguous hits to the same reference viral genome. If spacing of
# the hits on the bacterial genome roughly matches the spacing of these same regions on the viral genome, we merge the
# hits.
def assign_integration_ids(hit_list: List[QueryHit], overlap_tolerance_percent: float) -> Dict[int, List[QueryHit]]:
    #hit_list.sort(key=lambda x: (x.query_name, x.hit_name, find_position_for_strand_type(x.query_st, x.query_end, x.strand)))
    prev_hit = None
    # dict holding a list of hits for each individual integration, using integration ID as the key. Will be used later
    # to ask if an integration is considering full length based on whether
    integration_id_dict = {}
    integration_index = 0

    # for every hit in our list of detected viral sequences, compare to the previous hit in the list to determine if
    # they belong to the same integration
    for current_hit in hit_list:
        # if there is no previous hit, or if they aren't part of the same integration, iterate integration_index
        if prev_hit is None or not same_integration(prev_hit, current_hit, overlap_tolerance_percent):
            prev_hit = current_hit
            # we want the index to start at one, so it's fine to iterate before assigning
            integration_index += 1
            current_hit.integration_id = integration_index
            integration_id_dict[integration_index] = [current_hit]

        # otherwise, set to the same integration id as the previous hit
        else:
            current_hit.integration_id = integration_index
            integration_id_dict[integration_index].append(current_hit)

            prev_hit = current_hit

    return integration_id_dict


def set_integration_full_length(integration_dict: Dict[int, List[QueryHit]], full_threshold: float) -> None:
    # for each list of hits in an integration, if they cover enough of the reference viral genome to collectively
    # qualify as full length, set each hit to full_length = True
    for hit_list in integration_dict.values():
        if len(hit_list) > 1:
            percent_sum = 0.0
            for hit in hit_list:
                percent_sum += hit.get_percent_complete()

            if percent_sum >= full_threshold:
                for hit in hit_list:
                    hit.full_length = True


def sort_hit_list(hit_list: List[QueryHit]) -> None:
    # Sort the list, first by query sequence name (in case of multiple contigs) and then by start position relative to
    # the bacterial genome (so if the integration has - for strand, we want to use the end, which occurs first on the
    # bacterial genome). This will ensure that any hits that are part of the same integration are next to each other
    # in the list.
    hit_list.sort(key=lambda x: (x.query_name, x.hit_name, find_position_for_strand_type(x.query_st, x.query_end, x.strand)))


def overwrite_check(file_path: str, force: bool) -> None:
    if path.isfile(file_path) and not force:
        raise FileExistsError(
            f"Output file {file_path} already exists- either move or delete this file or enable --force")


def write_occurrence_json(json_file: TextIO, query_hits: List[QueryHit]) -> None:
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


def write_annotation_json(json_file: TextIO, query_hits: List[QueryHit], protein_annotations: Optional[Dict[str, str]],
                          occurrence_dict: Dict[str, List[int]], genome_path: str, verbose: bool) -> None:
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

    json_file.write(json.dumps(annotation_json, indent=INDENT_VAL))

def write_annotation_json_from_path(json_path: str, query_hits: List[QueryHit],
                                    protein_annotations: Optional[Dict[str, str]],
                                    occurrence_dict: Dict[str, List[int]], genome_path: str, force: bool,
                                    verbose: bool) -> None:
    overwrite_check(json_path, force)

    with open(json_path, "w") as json_file:
        write_annotation_json(json_file, query_hits, protein_annotations, occurrence_dict, genome_path, verbose)
        json_file.close()


def write_tsv(tsv: TextIO, seq_list: List[QueryHit]) -> None:
    # write header line first
    headers = ["Hit Name", "Description", "Accession", "E-Value", "Full Length",
               "Ref Seq Start", "Ref Seq End", "Ref Seq Length",
               "Bacterial Genome Name", "Query Sequence Name", "Match Start on Query Seq",
               "Match End on Query Seq", "Query Genome Length", "Strand", "Integration ID\n"]
    tsv.write("\t".join(headers))
    for viral_seq in seq_list:
        tsv.write(viral_seq.to_tsv_line())


def write_tsv_from_path(tsv_path: str, seq_list: List[QueryHit], force: bool) -> None:
    overwrite_check(tsv_path, force)

    with open(tsv_path, "w") as tsv:
        write_tsv(tsv, seq_list)


def set_hit_full_length(hit: QueryHit, threshold: float) -> None:
    percent_coverage_of_reference = hit.get_percent_complete()
    hit.full_length = hit.get_percent_complete() >= threshold


def parse_dfam_file(dfam_file: TextIO, genome_path: str, full_threshold: float, max_eval: float, verbose)\
        -> List[QueryHit]:
    hit_list = []
    # For the first hit, we don't know what the query genome length is. The QueryHit class contains a method that can
    # get the query genome length, but since we expect each run of table_parser.py to only deal with one query genome,
    # we can re-use the data one we get it via esl-seqstat
    genome_len = None
    for line_num, line in enumerate(dfam_file, 0):
        if line[0] == "#":
            pass
        else:
            line_list = line.split()

            # FraHMMER and nhmmscan disagree over what is the target and what is the query- we use FraHMMER's notation
            # TODO: Explanatory comments for each field
            hit_name = line_list[0]
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
                hit = QueryHit(hit_name, acc_id, query_name, evalue, ali_st, ali_en, genome_path, hmm_st, hmm_en,
                               hmm_len, strand, verbose, query_genome_len=genome_len)
                genome_len = hit.query_genome_len
                set_hit_full_length(hit, full_threshold)
                hit_list.append(hit)
            else:
                if verbose:
                    print(f"Excluding line {line_num}: e-value of {evalue} failed to pass maximum e-value threshold of "
                          f"{max_eval}")

    return hit_list


# TODO: Block comment
def parse_tbl_file(tbl_file: TextIO, genome_path: str, full_threshold: float, max_eval: float, verbose: bool,
                   annotations: Dict[str, str] = None) -> List[QueryHit]:
    hit_list = []
    genome_len = None
    for line_num, line in enumerate(tbl_file, 0):
        # skip comment lines starting with #
        if line[0] == "#":
            pass
        else:
            line_list = line.split()

            # FraHMMER and nhmmscan disagree over what is the target and what is the query- we use FraHMMER's notation
            # TODO: Explanatory comments for each field
            query_name = str(line_list[0])
            acc_id = str(line_list[1])
            hit_name = str(line_list[2])
            hmm_len = int(line_list[4])
            hmm_st = int(line_list[5]) # position of match relative to model
            hmm_en = int(line_list[6])
            ali_st = int(line_list[8])
            ali_en = int(line_list[9])
            evalue = float(line_list[12])
            description = "-"
            if annotations:
                description = annotations[hit_name]

            # TODO: comment on what strand is
            strand = "+"
            if hmm_en < hmm_st:
                strand = "-"

            if evalue <= max_eval:
                hit = QueryHit(hit_name, acc_id, query_name, evalue, ali_st, ali_en, genome_path, hmm_st, hmm_en,
                               hmm_len, strand, verbose, query_genome_len=genome_len, description=description)
                # all hits come from the same genome, so genome length is the same
                genome_len = hit.query_genome_len
                set_hit_full_length(hit, full_threshold)
                hit_list.append(hit)
            else:
                if verbose:
                    print(f"Excluding line {line_num}: e-value of {evalue} failed to pass maximum e-value threshold of "
                          f"{max_eval}")

    return hit_list


def parse_table_from_path(table_path: str, genome_path: str, full_threshold: float, max_eval: float,
                          table_mode: TABLE_MODE,  verbose: bool,
                          annotations: Dict[str, str] = None) -> List[QueryHit]:
    with open(table_path) as table_file:
        if verbose:
            print(f"Opening {table_path}...")

        if table_mode == "dfam":
            return parse_dfam_file(table_file, genome_path, full_threshold,  max_eval, verbose)
        elif table_mode == "tbl":
            return parse_tbl_file(table_file, genome_path, full_threshold, max_eval, verbose, annotations=annotations)
        else:
            raise ValueError("table_type must be either dfam or tbl")


def parse_protein_annotation_from_path(anno_tsv_path: str, verbose: bool) -> Dict[str, str]:
    if verbose:
        print(f"Opening {anno_tsv_path}...")

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


# TODO: subcommands, full length threshold %, minimum length for reporting, gap size %, overlap tolerance %
def parse_args(sys_args: List[str]) -> argparse.Namespace:
    default_eval = 1e-5
    parser = argparse.ArgumentParser(sys_args, description="Parses input .dfam or .tbl (from FraHMMER --tblout) file,"
                                    " extracting information about query hits detected in target genome")
    parser.add_argument("table_path", type=str, help="Path to input .dfam or .tbl file made up of query hits")
    parser.add_argument("genome_path", type=str, help="Path to target genome in .fasta format")
    parser.add_argument("output_tsv_path", type=str, help="Path to output .tsv file")
    parser.add_argument("output_json_path", type=str, help="Path to output .json file containing information on "
                                                           "nucleotide occurrence counts")
    parser.add_argument("table_type", type=str, choices=["dfam","tbl"], default="dfam", help="Which type of table is "
                        "being supplied as input, which must be dfam or tbl (default dfam)")
    parser.add_argument("json_type", type=str, choices=["integration","annotation"], default="integration", help="Sets"
                        " which type of .json will be produced: integration, which counts how many times a position in"
                        " a user-supplied viral genome has been detected in a specific bacterial genome, or annotation,"
                        " which produces .json describing a reference viral genome annotated with viral proteins."
                        " Default is integration")
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
    # TODO: Explanations
    args = parse_args(sys.argv[1:])
    table_path = args.table_path
    genome_path = args.genome_path
    tsv_path = args.output_tsv_path
    json_path = args.output_json_path
    table_mode = args.table_type
    json_mode = args.json_type
    occ_json_path = args.occurrence_json
    protein_annotation_path = args.annotation_tsv
    full_threshold = 0.7
    overlap_tolerance_percent = 0.01
    max_eval = args.max_evalue
    verbose = args.verbose
    force = args.force

    # check that inputs are legal
    if max_eval < 0:
        raise ValueError("--max_evalue must be used with an argument greater than or equal to 0")

    # integration mode
    if json_mode == "integration":
        # parse table for information on hits detected on query
        query_hits = parse_table_from_path(table_path, genome_path, full_threshold, max_eval, table_mode, verbose)
        # sort list to ensure that any hits from the same integration are next to each other
        sort_hit_list(query_hits)
        # examine sorted hits to determine if any of them are part of one integration broken up over multiple hits
        integration_id_dict = assign_integration_ids(query_hits, overlap_tolerance_percent)
        # check whether any integrations broken up over multiple hits cover enough of their reference viral genome to be
        # considered full length. If so, set each constituent hit to full_length = True
        set_integration_full_length(integration_id_dict, full_threshold)
        # write output
        write_tsv_from_path(tsv_path, query_hits, force)
        write_occurrence_json_from_path(json_path, query_hits, force)

    # annotation mode
    elif json_mode == "annotation":
        protein_annotations = None

        if protein_annotation_path:
            protein_annotations = parse_protein_annotation_from_path(protein_annotation_path, verbose)

        query_hits = parse_table_from_path(table_path, genome_path, full_threshold, max_eval, table_mode, verbose,
                                           protein_annotations)
        occurrence_dict = load_json(occ_json_path, verbose)
        write_tsv_from_path(tsv_path, query_hits, force)
        write_annotation_json_from_path(json_path, query_hits, protein_annotations, occurrence_dict, genome_path, force,
                                        verbose)


if __name__ == "__main__":
    _main()
