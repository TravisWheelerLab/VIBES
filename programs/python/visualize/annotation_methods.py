import re


class Genome:
    '''
    Describes a genome containing some matchs that we care about.

    Attributes:
    -----
    name: String
        Name of the genome.

    matches: List of Match objects
        A dictionary of matches to some sequence in the genome. Keys are the start position of the match in the genome, values are lists of Match objects
        (to cover the case of two matches having the same start position)
        '''
    def __init__(self, name):
        self.name = name
        self.matches = {}


class Match:
    '''
    Describes a match between some target sequence and some query sequence.

    Attributes:
    -----
    name: String
        Name of the target sequence.

    eval: float
        E-value of match, >= 0. Lower is better.

    hmm_st: int
        Where the match begins on the HMM. Because backwards matches occur, hmm_st is not necessarily smaller than hmm_en.

    hmm_en: int
        Where the match ends on the HMM. Because backwards matches occur, hmm_en is not necessarily larger than hmm_st.

    ali_st: int
        Where the match begins on the query genome. Because backwards matches occur, aliSt is not necessarily smaller than ali_en.

    ali_en: int
        Where the match ends on the query genome. Because backwards matches occur, ali_en is not necessarily larger than ali_st.

    acc_id: String
        Unique identifier of the target sequence.

    description: String
        Brief description of the target sequence.

    genome_length: int
        Length of the query genome.
    '''

    def __init__(self, name, eval, hmm_st, hmm_en, ali_st, ali_en, acc_id, description, genome_length=None):
        self.name = name
        self.eVal = eval
        self.hmmSt = hmm_st
        self.hmmEn = hmm_en
        self.aliSt = ali_st
        self.aliEn = ali_en
        self.accID = acc_id
        self.description = description
        self.genomeLength = genome_length


# .fasta format is a nightmare, but its header can roughly be divided into 2 parts: the first sequence of characters following >, up to the first space (the name)
# and everything afterward (the description). Since HMMER discards non-name parts of .fasta headers, we need to recover them from the original .fasta file to
# determine whether or not a sequnce corresponds to an integrase
def read_fasta_descs(fasta_path):
    desc_dict = {}

    with open(fasta_path, 'r') as fasta_file:
        for line in fasta_file:
            # if line starts with > when ignoring non-ASCII characters, capture everything before and after the first space character
            if re.match(r'(?a)>', line):
                regex_match = re.search(r'(?a)>(.+?) (.+?)\n', line)
                #if regex_match is None:
                    #print(fastaPath)
                    #print(line)

                name = regex_match.group(1)
                description = regex_match.group(2)

                desc_dict[name] = description

    return desc_dict


def detect_overlap(reg1_st, reg1_en, reg2_st, reg2_en):
    overlap = False

    # if reg1 starts before reg2 ends, and reg1 ends after reg2 starts, they overlap
    if reg1_st < reg2_en and reg1_en > reg2_st:
        overlap = True

    return overlap


def annotate_genome(prot_domtbl_dir, pfam_domtbl_dir, dfam_dir, prophage_name, min_eval, genome_length=None):
    '''
    Returns Genome object with its matches dictionary populated start location keys and lists of Match objects starting at that location. The contents of 
    matches are the annotation of the genome, ordered by a combination of start position relative to the genome and size (matches passing a length threshold
    occur before shorter, overlapping matches)
    '''

    # first, ensure that either protDomtblDir or pfamDomtblDir has been specified by user (this will allow us to provide genome length to Match objects generated
    # from .dfam files
    if prot_domtbl_dir is None and pfam_domtbl_dir is None:
        print("At least one of prot_domtbl_dir or pfam_domtbl_dir must be specified")
        exit()

    # if directory path string isn't empty, then read in information from file in that directory
    if prot_domtbl_dir:
        dom_tbl_path = "%s/%s.domtbl" % (prot_domtbl_dir, prophage_name)
        prot_domtbl_list = build_domtbl_list(dom_tbl_path, min_eval, "swissProt")

    if pfam_domtbl_dir:
        pfam_domtbl_path = "%s/%s.domtbl" % (pfam_domtbl_dir, prophage_name)
        pfam_domtbl_list = build_domtbl_list(pfam_domtbl_path, min_eval, "Pfam A")

    if dfam_dir:
        dfam_path = "%s/%s.dfam" % (dfam_dir, prophage_name)
        dfam_list = build_dfam_list(dfam_path, min_eval)
        if genome_length:
            set_dfam_genome_length(dfam_list, genome_length)
    else:
        dfam_list = []

    annotated_genome = Genome(prophage_name)
    # list of tuples containing protein start and end coordinates, relative to the viral genome
    prot_coord_list = []

    # currently, we expected .dfam format data to be recombinase and pseudogene hits. These are both (or were) protein-coding sequence, so we union
    # the list of .dfam Matchs with the list of protein .domtbl Matches
    for prot_match in (prot_domtbl_list + dfam_list):
        # to determine cases where protein domain annotations overlap with protein annotations in the prophage genomes,
        # we also set values corresponding to protein annotation coordinates to True in isOccupiedList.
        # Since isOccupiedList is 0-indexed, we subtract 1 from x_start
        x_start = prot_match.aliSt - 1
        x_end = prot_match.aliEn

        # record protein coordinate tuple in prot_coord_list. If end coordinate is smaller than start, flip them temporarily
        if x_end < x_start:
            temp_x_start = x_end
            temp_x_end = x_start
        else:
            temp_x_start = x_start
            temp_x_end = x_end

        prot_coord_list.append((temp_x_start, temp_x_end))

        if x_start in annotated_genome.matches:
            annotated_genome.matches[x_start].append(prot_match)
        else:
            value_list = [prot_match]
            annotated_genome.matches[x_start] = value_list

    for pfam_match in pfam_domtbl_list:
        x_start = pfam_match.aliSt - 1
        x_end = pfam_match.aliEn

        # should be False if no protein annotation lines overlap with domain annotation
        overlaps_with_protein = False

        for coord_tuple in prot_coord_list:
            # if x_end is less than x_start, temporarily flip them
            if x_end < x_start:
                temp_x_start = x_end
                temp_x_end = x_start
            else:
                temp_x_start = x_start
                temp_x_end = x_end

                if detect_overlap(temp_x_start, temp_x_end, coord_tuple[0], coord_tuple[1]):
                    overlaps_with_protein = True

        # We want to sort annotation lines by length to prioritize drawing longest lines closest to the x-axis, so we use line length as the key.
        # If key already in dictionary, append our list of line info to value (list of lists of line info)
        if not overlaps_with_protein:
            if x_start in annotated_genome.matches:
                annotated_genome.matches[x_start].append(pfam_match)
            else:
                value_list = [pfam_match]
                annotated_genome.matches[x_start] = value_list

    return annotated_genome


def set_dfam_genome_length(dfam_list, genome_length):
    for match in dfam_list:
        match.genome_length = genome_length


# read in information from .dfam format files
def build_dfam_list(dfam_path, min_eval):
    info_list = []

    # read in info from .dfam file
    with open(dfam_path, "r") as dfam_data:
        for line in dfam_data:
            # '#' char indicates a line doesn't contain data
            if line[0] != "#":
                # create a list to store this line's data
                data_list = line.split()
                join_string = " "

                # since split() gives us strings, we cast to the proper type
                match_name = data_list[0]
                i_evalue = float(data_list[4])
                hmm_from = int(data_list[6])
                hmm_to = int(data_list[7])
                ali_from = int(data_list[9])
                ali_to = int(data_list[10])
                description = join_string.join(data_list[14:])

                # temporary solution, since currently all expected .dfam files have had match names purposefully made short by me
                acc_id = match_name

                # we only want entries with e-value <= minimum (default 1e-5)
                if i_evalue <= min_eval:
                    match = Match(match_name, i_evalue, hmm_from, hmm_to, ali_from, ali_to, acc_id, description)
                    info_list.append(match)

    return info_list


# Read in the contents of a .domtbl file. Returns a list of lists of
# relevant data. Each list in the LoL corresponds to one line of the file
def build_domtbl_list(domtbl_path, min_eval, file_source):
    info_list = []

    # read in domTblPath info as read-only
    with open(domtbl_path, "r") as domtbl_data:
        for line in domtbl_data:
            # '#' char indicates a line doesn't contain data
            if line[0] != "#":
                # create a list to store this line's datallinux terminal is something a fileinux terminal is something a file
                data_list = line.split()
                join_string = " "

                # since split() gives us strings, we cast to the proper type
                domain_name = data_list[0]
                tlen = int(data_list[2])
                i_evalue = float(data_list[12])  # i-evalue is domain-specific evalue
                hmm_from = int(data_list[15])
                hmm_to = int(data_list[16])
                ali_from = int(data_list[19])
                ali_to = int(data_list[20])
                description = join_string.join(data_list[27:])

                if file_source == "swissProt":
                    # use regex to extract accession ID from domain_name
                    name_search = re.search(r'\|(.+?)\|', domain_name)
                    acc_id = name_search[1]
                else:
                    acc_id = data_list[1]

                # we only want entries with e-value <= minimum (default 1e-5)
                if i_evalue <= min_eval:
                    match = Match(domain_name, i_evalue, hmm_from, hmm_to, ali_from, ali_to, acc_id, description, tlen)
                    info_list.append(match)

    return info_list
