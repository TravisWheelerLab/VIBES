import matplotlib.pyplot as plt
import numpy as np
import argparse
import re
from os import walk, path
import sys
from matplotlib import lines, cm
from matplotlib.colors import ListedColormap
from generate_domtbls import DomTblGen


# read in all files in total counts dir via walk
# for each file, populate an array in which each index corresponds to a line
# plot as line graph
def plotList(countList, prophageName, outputDir, domTblDir, hmmStatDict):

    print(prophageName)
    # Y depth constant. Used to fiddle with how far below the x-axis lines are drawn
    YCONST = -.1
    STACKCONST = .1

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    ax.plot(countList)

    domTblPath = "%s%s.domtbl" % (domTblDir, prophageName)
    domTblList = buildDomTblList(domTblPath)

    '''if (len(countList) == domTblList[0][1]):
        print("Lengths match!")



    print(prophageName) # TODO investigate why len(countList) and tlen disagree
    print(countList)
    print("Countlist: %d" % len(countList))
    print("Tlen: %s" % domTblList[0][1])
    exit()'''

    # CREATE LINE DEPTH STUFF SOMEWHERE AROUND HERE
    # In case protein domains overlap, we want to stagger the drawing of
    # the lines that represent them. To accomplish this, we create a list
    # of booleans set to false. When we insert a line on a "layer," we set
    # all values along that line on its layer to true. To determine which
    # layer we should draw a line on, we start at the lowest layer, and draw
    # on the first layer with enough free space for the entire line plus some
    # amount of buffer

    # create list of lists filled with False
    depthList = []
    for i in domTblList:
        depthList.append([False] * len(nucleotideList))

    # find max height of chart. Useful for determining y-values below
    YMAX = max(countList)
    domainLines = []
    descStrings = []

    for lineList in domTblList:
        # grab x and y data, from alicoords and line depth stuff respectively.
        # subtract 1 from x-coords to adjust for 0-indexed graph
        xStart = lineList[5] - 1
        xEnd = lineList[6] - 1
        description = lineList[-1]

        lineLayer = 0
        while (True in depthList[lineLayer][xStart:(xEnd + 1)]):
            lineLayer += 1

        #print(lineLayer)

        index = xStart

        while (index < xEnd + 1):
            depthList[lineLayer][index] = True
            #print("Index: %d" % index)
            index += 1

        #print(depthList[lineLayer])

        print("done")

        y = (YMAX * YCONST) - ((STACKCONST * YMAX) * lineLayer)

        line = lines.Line2D(np.array([xStart, xEnd]), np.array([y, y]), clip_on=False, linewidth=3)

        ax.add_line(line)

        # create annotation for line, based on line's center point and offset to be under line
        xCenter = xStart + abs((xEnd - xStart) / 2)
        ax.annotate(description, xy=(xCenter, y - (y * STACKCONST * 2)), annotation_clip=False, horizontalalignment='center', fontsize=10)

        domainLines.append(line)
        descStrings.append(description)

    # plot line. Use number of domain lines to determine how many colors we need
    # set clip to false
    # plot

    # Create colormap. Credit to https://stackoverflow.com/questions/4971269/how-to-pick-a-new-color-for-each-plotted-line-within-a-figure-in-matplotlib
    colors = iter(cm.jet(np.linspace(0, 1, len(domainLines))))
    for i in range(len(domainLines)):
        currentColor = next(colors)
        domainLines[i].set_color(currentColor)

    ax.set(xlabel='Nucleotide Index', ylabel='Number Found', title=prophageName)
    # credit to https://stackoverflow.com/questions/7125009/how-to-change-legend-size-with-matplotlib-pyplot for legend settings
    #ax.legend(domainLines, descStrings, loc=1, prop={'size': 6})
    #ax.grid()

    plt.tight_layout()
    fig.savefig("%s/%s.png" % (outputDir,prophageName), transparent=True)
    plt.show()
    plt.close()


# Read in the contents of a .domtbl file. Returns a list of lists of
# relevant data. Each list in the LoL corresponds to one line of the file
def buildDomTblList(domTblPath):
    infoList = []

    # read in domTblPath info as read-only
    with open(domTblPath, "r") as domTblData:
        for line in domTblData:
            # '#' char indicates a line doesn't contain data
            if(line[0] != "#"):
                # create a list to store this line's data
                lineList = []
                dataList = line.split()

                # since we split domTblData along whitespace, we need to reconstruct
                # the description string at the end of the table line. This string
                # always begins at index 27 and runs until the end of the list
                descString = " ".join(dataList[27:])

                # since split() gives us strings, we cast to the proper type
                domainName = dataList[0]
                tlen = int(dataList[2])
                iEvalue = float(dataList[12])  # i-Evalue is domain-specific Evalue
                hmmFrom = int(dataList[15])
                hmmTo = int(dataList[16])
                aliFrom = int(dataList[19])
                aliTo = int(dataList[20])

                # we only want entries with evalue 10^-5 or better
                if (iEvalue <= 1e-5):
                    print("Passing eVal: %s" % iEvalue)
                    lineList.append(domainName)
                    lineList.append(tlen)
                    lineList.append(iEvalue)
                    lineList.append(hmmFrom)
                    lineList.append(hmmTo)
                    lineList.append(aliFrom)
                    lineList.append(aliTo)
                    lineList.append(descString)

                    infoList.append(lineList)

    return infoList


# reads in model length field ('M') for each domain in a .hmmstattbl file
# returns a dictionary containing all domains in file, with domain name key
# and model length value
def buildHmmStatDict(hmmStatPath):
    statDict = {}

    with open(hmmStatPath, "r") as hmmStatData:
        for line in hmmStatData:
            # '#' char indicates a line doesn't contain data, 'not line'
            # checks for the empty line that occurs at EOF
            if(line[0] != "#" and not line):
                dataList = line.split()
                print("Line: %s" % (line))
                domName = dataList[1]
                modelLength = dataList[6]

                statDict[domName] = modelLength

    return statDict


# credit to Viktor Kerkez in: https://stackoverflow.com/questions/18160078/how-do-you-write-tests-for-the-argparse-portion-of-a-python-module
def parseArgs(sysArgs):
    parser = argparse.ArgumentParser(sysArgs)
    parser.add_argument("countDir", help="Path to directory of nucleotide count files. Expected format is [prophageName]Chart.txt")
    parser.add_argument("domTblDir", help="Path to directory of .domtbl files. Expected format is [prophageName].domtbl")
    parser.add_argument("hmmStat", help="Path to a hmmstat .tbl file")
    parser.add_argument("outputDir", help="Path to directory to store output .png files")
    parser.add_argument("--force", help="If output files already exist, overwrite them", action="store_true")

    return parser.parse_args()


if __name__ == "__main__":
    args = parseArgs(sys.argv[1:])
    countDir = args.countDir
    hmmStatPath = args.hmmStat

    # builds dictionary containing model length for each domain
    hmmStatDict = buildHmmStatDict(hmmStatPath)

    # credit to pycruft in https://stackoverflow.com/questions/3207219/how-do-i-list-all-files-of-a-directory for code for grabbing file paths
    filePaths = []
    for (dirpath, dirnames, files) in walk(countDir):
        filePaths.extend(files)

    filePaths.sort()

    # credit to Tim Anderson for advice on how to open files:
    # runs through all paths in the directory, opening as read only
    for filePath in filePaths:
        nucleotideList = []

        with open("%s/%s" % (countDir, filePath), "r") as fileData:
            # skip header line
            fileData.readline()
            for line in fileData:
                line.rstrip("\n")
                nucleotideList.append(int(line))

        regMatch = re.match(r'(.+?)Chart\.txt', filePath)
        # extract name from file path
        prophageName = regMatch.group(1)

        #print(len(nucleotideList))

        plotList(nucleotideList, prophageName, args.outputDir, args.domTblDir, hmmStatDict)