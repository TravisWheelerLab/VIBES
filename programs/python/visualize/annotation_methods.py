import re


def genAnnotationDict(protDomtblDir, pfamDomtblDir, dfamDir, prophageName, minEval):
    genomeLength = None

    # first, ensure that either protDomtblDir or pfamDomtblDir points to a directory (this will allow us to find genome length)
    if protDomtblDir is None and pfamDomtblDir is None:
        print("At least one of protDomtblDir or pfamDomtblDir must be specified")
        exit()

    # if directory path string isn't empty, then read in information from file in that directory
    if protDomtblDir:
        domTblPath = "%s/%s.domtbl" % (protDomtblDir, prophageName)
        protDomtblLists = buildDomtblList(domTblPath, minEval, "swissProt")

        # we grab genomeLength from first entry of list of line info lists for use in buildDfamList
        genomeLength = protDomtblLists[0][1]

    if pfamDomtblDir:
        pfamDomtblPath = "%s/%s.domtbl" % (pfamDomtblDir, prophageName)
        pfamDomtblLists = buildDomtblList(pfamDomtblPath, minEval, "Pfam A")

        # we grab genomeLength from first entry of list of line info lists for use in buildDfamList
        genomeLength = pfamDomtblLists[0][1]

    if dfamDir:
        dfamPath = "%s/%s.dfam" % (dfamDir, prophageName)
        dfamLists = buildDfamList(dfamPath, minEval, genomeLength)

    annotationDict = {}
    # list of booleans with a length equivalent to genome being plotted. We want to avoid annotating these images with protein domains
    # in locations already annotated by proteins, so indexes annotated by protein data will be set to true. If a domain's annotation overlaps
    # with any of these true values, we don't include it
    isOccupiedList = [False] * genomeLength

    # currently, we expected .dfam format data to be recombinase and pseudogene matches
    for protLineList in (protDomtblLists + dfamLists):
        # to determine cases where protein domain annotations overlap with protein annotations
        # we also set values corresponding to protein annotation coordinates to True in isOccupiedList.
        # Since isOccupiedList is 0-indexed, we subtract 1 from xStart
        xStart = protLineList[5] - 1
        xEnd = protLineList[6]

        for index in range(xStart, xEnd):
            isOccupiedList[index] = True

        if xStart in annotationDict:
            annotationDict[xStart].append(protLineList)
        else:
            valueList = [protLineList]
            annotationDict[xStart] = valueList

    for pfamLineList in pfamDomtblLists:
        xStart = pfamLineList[5] - 1
        xEnd = pfamLineList[6]

        # should be False if no protein annotation lines overlap with domain annotation
        overlapsWithProtein = False

        for index in range(xStart, xEnd):
            if isOccupiedList[index]:
                overlapsWithProtein = True

        # We want to sort annotation lines by length to prioritize drawing longest lines closest to the x-axis, so we use line length as the key.
        # If key already in dictionary, append our list of line info to value (list of lists of line info)
        if not overlapsWithProtein:
            if xStart in annotationDict:
                annotationDict[xStart].append(pfamLineList)
            else:
                valueList = []
                valueList.append(pfamLineList)
                annotationDict[xStart] = valueList

    return annotationDict


# read in information from .dfam format files
def buildDfamList(dfamPath, minEval, genomeLength):
    infoList = []

    # read in info from .dfam file
    with open(dfamPath, "r") as dfamData:
        for line in dfamData:
            # '#' char indicates a line doesn't contain data
            if(line[0] != "#"):
                # create a list to store this line's data
                lineList = []
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
                    lineList.append(matchName)
                    lineList.append(genomeLength)
                    lineList.append(iEvalue)
                    lineList.append(hmmFrom)
                    lineList.append(hmmTo)
                    lineList.append(aliFrom)
                    lineList.append(aliTo)
                    lineList.append(accID)
                    lineList.append(description)

                    infoList.append(lineList)

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
                # create a list to store this line's data
                lineList = []
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
                    lineList.append(domainName)
                    lineList.append(tlen)
                    lineList.append(iEvalue)
                    lineList.append(hmmFrom)
                    lineList.append(hmmTo)
                    lineList.append(aliFrom)
                    lineList.append(aliTo)
                    lineList.append(accID)
                    lineList.append(description)

                    infoList.append(lineList)

    return infoList
