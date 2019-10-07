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

    eVal: float
        E-value of match, >= 0. Lower is better.

    hmmSt: int
        Where the match begins on the HMM. Because backwards matches occur, hmmSt is not necessarily smaller than hmmEn.

    hmmEn: int
        Where the match ends on the HMM. Because backwards matches occur, hmmSt is not necessarily smaller than hmmEn.

    aliSt: int
        Where the match begins on the query genome. Because backwards matches occur, aliSt is not necessarily smaller than aliEn.

    accID: String
        Unique identifier of the target sequence.

    description: String
        Brief description of the target sequence.

    genomeLength: int
        Length of the query genome.
    '''

    def __init__(self, name, eVal, hmmSt, hmmEn, aliSt, aliEn, accID, description, genomeLength=None):
        self.name = name
        self.eVal = eVal
        self.hmmSt = hmmSt
        self.hmmEn = hmmEn
        self.aliSt = aliSt
        self.aliEn = aliEn
        self.accID = accID
        self.description = description
        self.genomeLength = genomeLength


def detectOverlap(reg1St, reg1En, reg2St, reg2En):
    overlap = False

    # if reg1 starts before reg2 ends, and reg1 ends after reg2 starts, they overlap
    if reg1St < reg2En and reg1En > reg2St:
        overlap = True

    return overlap


def annotateGenome(protDomtblDir, pfamDomtblDir, dfamDir, prophageName, minEval, genomeLength=None):
    '''
    Returns Genome object with its matches dictionary populated start location keys and lists of Match objects starting at that location. The contents of 
    matches are the annotation of the genome, ordered by a combination of start position relative to the genome and size (matches passing a length threshold
    occur before shorter, overlapping matches)
    '''

    # first, ensure that either protDomtblDir or pfamDomtblDir has been specified by user (this will allow us to provide genome length to Match objects generated
    # from .dfam files
    if protDomtblDir is None and pfamDomtblDir is None:
        print("At least one of protDomtblDir or pfamDomtblDir must be specified")
        exit()

    # if directory path string isn't empty, then read in information from file in that directory
    if protDomtblDir:
        domTblPath = "%s/%s.domtbl" % (protDomtblDir, prophageName)
        protDomtblList = buildDomtblList(domTblPath, minEval, "swissProt")

    if pfamDomtblDir:
        pfamDomtblPath = "%s/%s.domtbl" % (pfamDomtblDir, prophageName)
        pfamDomtblList = buildDomtblList(pfamDomtblPath, minEval, "Pfam A")

    if dfamDir:
        dfamPath = "%s/%s.dfam" % (dfamDir, prophageName)
        dfamList = buildDfamList(dfamPath, minEval)
        if genomeLength:
            setDfamGenomeLength(dfamList, genomeLength)

    annotatedGenome = Genome(prophageName)
    # list of tuples containing protein start and end coordinates, relative to the viral genome
    protCoordList = []

    # currently, we expected .dfam format data to be recombinase and pseudogene hits. These are both (or were) protein-coding sequence, so we union
    # the list of .dfam Matchs with the list of protein .domtbl Matches
    for protMatch in (protDomtblList + dfamList):
        # to determine cases where protein domain annotations overlap with protein annotations in the prophage genomes,
        # we also set values corresponding to protein annotation coordinates to True in isOccupiedList.
        # Since isOccupiedList is 0-indexed, we subtract 1 from xStart
        xStart = protMatch.aliSt - 1
        xEnd = protMatch.aliEn

        # record protein coordinate tuple in protCoordList. If end coordinate is smaller than start, flip them temporarily
        if xEnd < xStart:
            tempXStart = xEnd
            tempXEnd = xStart
        else:
            tempXStart = xStart
            tempXEnd = xEnd

        protCoordList.append((tempXStart, tempXEnd))

        if xStart in annotatedGenome.matches:
            annotatedGenome.matches[xStart].append(protMatch)
        else:
            valueList = [protMatch]
            annotatedGenome.matches[xStart] = valueList

    for pfamMatch in pfamDomtblList:
        xStart = pfamMatch.aliSt - 1
        xEnd = pfamMatch.aliEn

        # should be False if no protein annotation lines overlap with domain annotation
        overlapsWithProtein = False

        for coordTuple in protCoordList:
            # if xEnd is less than xStart, temporarily flip them
            if xEnd < xStart:
                tempXStart = xEnd
                tempXEnd = xStart
            else:
                tempXStart = xStart
                tempXEnd = xEnd

                if detectOverlap(tempXStart, tempXEnd, coordTuple[0], coordTuple[1]):
                    overlapsWithProtein = True

        # We want to sort annotation lines by length to prioritize drawing longest lines closest to the x-axis, so we use line length as the key.
        # If key already in dictionary, append our list of line info to value (list of lists of line info)
        if not overlapsWithProtein:
            if xStart in annotatedGenome.matches:
                annotatedGenome.matches[xStart].append(pfamMatch)
            else:
                valueList = [pfamMatch]
                annotatedGenome.matches[xStart] = valueList

    return annotatedGenome

def setDfamGenomeLength(dfamList, genomeLength):
    for match in dfamList:
        match.genomeLength = genomeLength


# read in information from .dfam format files
def buildDfamList(dfamPath, minEval):
    infoList = []

    # read in info from .dfam file
    with open(dfamPath, "r") as dfamData:
        for line in dfamData:
            # '#' char indicates a line doesn't contain data
            if(line[0] != "#"):
                # create a list to store this line's data
                dataList = line.split()
                joinString = " "

                # since split() gives us strings, we cast to the proper type
                matchName = dataList[0]
                iEvalue = float(dataList[4])
                hmmFrom = int(dataList[6])
                hmmTo = int(dataList[7])
                aliFrom = int(dataList[9])
                aliTo = int(dataList[10])
                description = joinString.join(dataList[14:])

                # temporary solution, since currently all expected .dfam files have had match names purposefully made short by me
                accID = matchName

                # we only want entries with e-value <= minimum (default 1e-5)
                if (iEvalue <= minEval):
                    match = Match(matchName, iEvalue, hmmFrom, hmmTo, aliFrom, aliTo, accID, description)
                    infoList.append(match)

    return infoList


# Read in the contents of a .domtbl file. Returns a list of lists of
# relevant data. Each list in the LoL corresponds to one line of the file
def buildDomtblList(domTblPath, minEval, fileSource):
    infoList = []

    # read in domTblPath info as read-only
    with open(domTblPath, "r") as domTblData:
        for line in domTblData:
            # '#' char indicates a line doesn't contain data
            if(line[0] != "#"):
                # create a list to store this line's datallinux terminal is something a fileinux terminal is something a file
                dataList = line.split()
                joinString = " "

                # since split() gives us strings, we cast to the proper type
                domainName = dataList[0]
                tlen = int(dataList[2])
                iEvalue = float(dataList[12])  # i-Evalue is domain-specific Evalue
                hmmFrom = int(dataList[15])
                hmmTo = int(dataList[16])
                aliFrom = int(dataList[19])
                aliTo = int(dataList[20])
                description = joinString.join(dataList[27:])

                if fileSource == "swissProt":
                    # use regex to extract accession ID from domainName
                    nameSearch = re.search(r'\|(.+?)\|', domainName)
                    accID = nameSearch[1]
                else:
                    accID = dataList[1]

                # we only want entries with e-value <= minimum (default 1e-5)
                if (iEvalue <= minEval):
                    match = Match(domainName, iEvalue, hmmFrom, hmmTo, aliFrom, aliTo, accID, description, tlen)
                    infoList.append(match)

    return infoList
