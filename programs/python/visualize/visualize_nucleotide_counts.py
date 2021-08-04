import matplotlib.pyplot as plt
import numpy as np
import argparse
import re
from os import walk
from annotation_methods import annotateGenome, Genome
import sys
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import lines, cm
import json

# creates json file with necessary information to draw graph in SODA rather than matplotlib. This .json
# will contain the prophage name, length, index occurence counts, and protein annotation information
def to_json(countList, prophageName, jsonDir, protDomtblDir, pfamDomtblDir, dfamDir, minEval):
    INDENT_VALUE = 4
    genomeLength = len(countList)
    # returns Genome object, containing dictionary populated with Match objects
    annotationGenome = annotateGenome(protDomtblDir, pfamDomtblDir, dfamDir, prophageName, minEval, genomeLength)

    # begin populating json entries
    jsonDict = {"prophageName": prophageName,
                "occurrences": countList}

    # list that will hold protein annotation information
    protAnnoList = []
    # the starting index serves as the key for the matches dict, which has a
    # value made up of a list of Match objects whose annotation starts at that
    # index
    for startIndex in annotationGenome.matches:
        for match in annotationGenome.matches[startIndex]:
            matchDict =    {"name": match.name,
                            "ID": match.accID, 
                            "genome_st": match.aliSt,
                            "genome_en": match.aliEn,
                            "gene_st" : match.hmmSt,
                            "gene_en" : match.hmmEn,
                            "gene_len" : match.geneLength,
                            "evalue": match.eVal,
                            "desc": match.description}

            protAnnoList.append(matchDict)

    jsonDict["annotations"] = protAnnoList

    jsonFilePath = "%s/%s.json" % (jsonDir, prophageName)

    with open(jsonFilePath, 'w') as jsonFile:
        jsonFile.write(json.dumps(jsonDict, indent=INDENT_VALUE))
        jsonFile.close()


# read in all files in total counts dir via walk
# for each file, populate an array in which each index corresponds to a line
# plot as line graph
def drawPlot(countList, prophageName, outputDir, protDomtblDir, pfamDomtblDir, dfamDir, minEval):
    # Minimum length a line should be to be added to 'long-line' dictionary
    LONG_ANNO_CONST = .04
    genomeLength = len(countList)
    longAnnoLength = LONG_ANNO_CONST * genomeLength
    annotationGenome = Genome(prophageName)

    print(prophageName)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax.set(xlabel='Position', ylabel='Occurrences', title=prophageName)
    ax.plot(countList)

    annotationGenome = annotateGenome(protDomtblDir, pfamDomtblDir, dfamDir, prophageName, minEval, genomeLength)
    longAnnoDict = {}
    shortAnnoDict = {}

    # Below, we generate our long annotation list and short annotation list. We create two lists because we want longer lines to be
    # drawn closer to the x-axis, but we also want lines drawn in order of start position as a secondary priority. To accomplish this,
    # we first plot lines we consider long enough to be especially interesting in order of start pos, then shorter ones in order of start pos
    for key in annotationGenome.matches:
        for match in annotationGenome.matches[key]:
            annoLineLength = abs(match.aliSt - match.aliEn) + 1

            # if long enough to be considered a long line, we add to long line dictionary (for which each key is the start position and each value is a list of line lists)
            if(annoLineLength >= longAnnoLength):
                if key in longAnnoDict:
                    longAnnoDict[key].append(match)
                else:
                    valueList = []
                    valueList.append(match)
                    longAnnoDict[key] = valueList
            else:
                if key in shortAnnoDict:
                    shortAnnoDict[key].append(match)
                else:
                    valueList = []
                    valueList.append(match)
                    shortAnnoDict[key] = valueList


    # find max height of chart. Used in plotting lines and labels
    YMAX = max(countList)

    # Generate an array of colors, create a counter to loop through them
    colors = cm.Dark2(np.linspace(0, 1, 8))
    colorCount = 0

    # sort keys of dictionary containing lists of line information. Keys correspond to starting index. Put all keys into one list and sort to make coloring easier
    sortedLongKeys = sorted(longAnnoDict)
    sortedShortKeys = sorted(shortAnnoDict)

    # depthList is intended to determine when annotations overlap and must be vertically staggered. To accomplish this, we populate depthList
    # with inner lists, which are each the length of the x-axis and filled with False. There should be one inner list per annotation, to account
    # for cases where all annotation lines overlap
    depthList = []

    # we store all plotted line2D objects in a list, so we can assign color later based on starting position on x-axis
    line2DList = []

    # draw lines
    plotAnnotationDict(fig, ax, longAnnoDict, sortedLongKeys, depthList, YMAX, line2DList)
    plotAnnotationDict(fig, ax, shortAnnoDict, sortedShortKeys, depthList, YMAX, line2DList)

    # sort lines by xStart
    line2DList = sorted(line2DList, key=lambda line: line.get_xdata()[0])

    for line2D in line2DList:
        # assign color to line based on colorCount. If colorCount exceeds size of colors list, reset to 0
        color = colors[colorCount]
        line2D.set_color(color)
        colorCount += 1

        if (colorCount == len(colors)):
            colorCount = 0

    plt.tight_layout()
    plt.xlim(0, genomeLength)
    plt.ylim(bottom=0)
    pp = PdfPages("%s/%s.pdf" % (outputDir, prophageName))
    pp.savefig(bbox_inches='tight')
    pp.close()
    #plt.show()
    plt.close()


# Accepts a MatPlotLib axis, a dictionary, and a list of dictionary keys. Each value in the dictionary is expected to be a list of lists,
# where each inner list contains information necessary to draw and label each line. The order of keys in keyList determines the order in which
# lines are drawn
def plotAnnotationDict(fig, ax, lineListDict, keyList, depthList, YMAX, line2DList):
    # How far below the x-axis we want to start drawing lines and labels, in pixels
    Y_CONST = 100  # how far below x-axis in pixels to begin drawing lines
    # Offset between lines on different vertical layers (depth layers)
    STACK_CONST = 30  # offset between lines in display units (pixels?)
    # How far above annotation lines labels are placed
    LABEL_OFFSET_CONST = .014
    # How much space each character in a label is estimated to take

    for key in keyList:
        valueList = lineListDict[key]

        for match in valueList:
            # grab x and y data from alicoords and line depth stuff respectively.
            # subtract 1 from x-coords to adjust for 0-indexed graph
            xStart = match.aliSt - 1
            xEnd = match.aliEn - 1
            accID = match.accID
            genomeLength = match.genomeLength

            # To determine width of annotation text label, we place the label on the plot and draw it, generating its size. We then get its coordinates and
            # convert them into data units rather than display units. We use this to determine whether or not a line is long enough to prevent label overlap,
            # or if we should pad out its ends (padding only affects line and label placement, and is not represented on the plot itself). Credit to
            # https://stackoverflow.com/questions/41267733/getting-the-coordinates-of-a-matplotlib-annotation-label-in-figure-coordinates
            ann = ax.annotate(accID, xy=(0, 0), annotation_clip=False, horizontalalignment='center', fontsize=10)

            fig.canvas.draw()

            dataInv = ax.transData.inverted()  # convert from display units (pixels?) to data units

            # cerdit to jklymak at https://stackoverflow.com/questions/18405652/savefig-pgf-runtimeerror-cannot-get-window-extent-w-o-renderer for
            # workaround to 'Cannot get window extent w/o renderer' error
            box = mpl.text.Text.get_window_extent(ann, renderer=fig.canvas.get_renderer())
            textAxCoords = dataInv.transform(box.extents)
            textWidth = textAxCoords[2] - textAxCoords[0]

            # CREATE LINE DEPTH STUFF SOMEWHERE AROUND HERE
            # In case protein domains overlap, we want to stagger the drawing of
            # the lines that represent them. To accomplish this, we create a list
            # of booleans set to false. When we insert a line on a "layer," we set
            # all values along that line on its layer to true. To determine which
            # layer we should draw a line on, we start at the lowest layer, and draw
            # on the first layer with enough free space for the entire line plus some
            # amount of buffer

            depthList.append([False] * genomeLength)

            # Sometimes, a sequence's start location occurs later than its end location. In this case, we swap the two to improve readability
            if (xEnd < xStart):
                temp = xStart
                xStart = xEnd
                xEnd = temp

            lineLength = xEnd - xStart + 1
            # Some annotation lines are long enough that we don't need to pad out their ends to prevent label overlap. In these cases,
            # skip adding padding
            if(lineLength < textWidth):
                # determine size of buffer to append to ends of protein x-axis coordinates. This should
                # improve protein line placement and graph readability by staggering lines that
                # occur near to each other, but don't overlap
                padding = int((textWidth - lineLength) / 2)
                paddedStart = xStart - padding
                paddedEnd = xEnd + padding

                # if adding padding extends end values to indexes that don't exist, trim back to 0 or last valid index
                if (paddedStart < 0):
                    paddedStart = 0

                if (paddedEnd > (genomeLength - 1)):
                    paddedEnd = genomeLength - 1
            else:
                paddedStart = xStart
                paddedEnd = xEnd

            # check if layer, starting with first, is occupied. If yes, increment layer value (corresponds to moving down on the graph,
            # preventing overlaps)
            lineLayer = 0
            index = paddedStart

            while (index <= paddedEnd):
                if depthList[lineLayer][index] is True:
                    # if a part of layer is occupied, reset progress to ensure beginning of next layer is free
                    index = paddedStart
                    lineLayer += 1
                index += 1

            #print(lineLayer)

            index = paddedStart

            while (index <= paddedEnd):
                depthList[lineLayer][index] = True
                #print("Index: %d" % index)
                index += 1

            #print(depthList[lineLayer])

            #yOffset = dataInv.transform(yOffset)

            #axesTrans = ax.transAxes  # convert from axes units (0-1) to display units

            # yOffset = axesTrans.transform((0, -.1725))
            xAxisPos = ax.transData.transform([0, 0])  # we use this to determine location of origin, which has height of x-axis, in pixels

            yPosList = [0, xAxisPos[1] - Y_CONST - (lineLayer * STACK_CONST)]

            yPosArray = dataInv.transform(yPosList)  # transform back to data units

            y = yPosArray[1]

            #y = (firstLayerLoc[1] + (lineLayer * StackConstUnits[1]))

            line = lines.Line2D(np.array([xStart, xEnd]), np.array([y, y]), clip_on=False, linewidth=3)
            ax.add_line(line)
            line2DList.append(line)

            # determine where we want label and set its position to new location
            xCenter = xStart + abs((xEnd - xStart) / 2)

            labelOffset = YMAX * LABEL_OFFSET_CONST
            annCoords = (xCenter, y + labelOffset)
            ann.set_position(annCoords)

            # To cut down on guesswork and ensure we don't have text annotation labels overlapping regardless of plot dimensions, we get the dimensions of
            # the text itself from the plot. To do this, we must draw the plot and then get the coordinates of each of its corners. We then convert the coordinates,
            # which are in display units, into figure axis units.
            # credit to https://stackoverflow.com/questions/41267733/getting-the-coordinates-of-a-matplotlib-annotation-label-in-figure-coordinates


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
                #print("Line: %s" % (line))
                domName = dataList[1]
                modelLength = dataList[6]

                statDict[domName] = modelLength

    return statDict


# credit to Viktor Kerkez in: https://stackoverflow.com/questions/18160078/how-do-you-write-tests-for-the-argparse-portion-of-a-python-module
def parseArgs(sysArgs):
    parser = argparse.ArgumentParser(sysArgs)
    parser.add_argument("count_dir", help="Path to directory of nucleotide count files. Expected format is [prophageName]Chart.txt")
    parser.add_argument("output_dir", help="Path to directory to store output .png or .json files")
    parser.add_argument("--swissprot_domtbl", help="Path to directory of .domtbl files resulting from a search against a database of Swissprot proteins. Expects Swissprot entries to begin with '>sp|[AccID]|', and for file names to be [prophageName].domtbl.")
    parser.add_argument("--pfam_domtbl", help="Path to directory of .domtbl files resulting from a search against a Pfam A database. AccIDs should be in the 'accession' column, and file names are expected to be [prophageName].domtbl. ")
    parser.add_argument("--dfam_dir", help="Path to directory of .dfam files annotating viral genomes. Expects .dfam files to contain matches of protein DNA sequences to viral genomes. Expects file names to be [prophageName].dfam")
    parser.add_argument("--force", help="If output files already exist, overwrite them", action="store_true")
    parser.add_argument("--json", help="Output .json files with necessary information to create plots, rather than .pngs", action="store_true")
    parser.add_argument("--min_eval", type=float, default=1e-5, help="Minimum e-value required for a .domtbl entry to be included. Default is 1e-5")

    return parser.parse_args()


if __name__ == "__main__":
    args = parseArgs(sys.argv[1:])
    countDir = args.count_dir
    minEval = args.min_eval
    protDomtblDir = args.swissprot_domtbl
    pfamDomtblDir = args.pfam_domtbl
    dfamDir = args.dfam_dir
    outputDir = args.output_dir
    jsonBool = args.json

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

        '''print(prophageName)
        print(len(nucleotideList))
        print("Hey, fix that you commented out plotList()\n")'''

        if jsonBool:
            to_json(nucleotideList, prophageName, outputDir, protDomtblDir, pfamDomtblDir, dfamDir, minEval)
        else:
            drawPlot(nucleotideList, prophageName, outputDir, protDomtblDir, pfamDomtblDir, dfamDir, minEval)
