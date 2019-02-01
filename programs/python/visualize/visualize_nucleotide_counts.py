import matplotlib.pyplot as plt
import numpy as np
import argparse
import re
from os import walk
from matplotlib import collections as col

# read in all files in total counts dir via walk
# for each file, populate an array in which each index corresponds to a line
# plot as line graph


def plotList(list, prophageName, outputDir):
    # Credit to https://matplotlib.org/gallery/lines_bars_and_markers/simple_plot.html#sphx-glr-gallery-lines-bars-and-markers-simple-plot-py
    # for structure of plotting code
    fig, ax = plt.subplots()
    ax.plot(list)

    maximum = max(list)
    yValue = maximum * -.13
    length = len(list)
    xValue = length/3

    lines = [[(0, yValue), (xValue, yValue)]]
    lc = col.LineCollection(lines)
    lc.set_clip_on(False)
    ax.add_collection(lc)

    ax.set(xlabel='Nucleotide Index', ylabel='Number Found', title=prophageName)
    ax.grid()

    fig.savefig("%s/%s.png" % (outputDir,prophageName))
    #plt.show()
    plt.close()

    def readInfo(domTblPath, hmmStatPath):
        infoDict = {}
        # read in domTblPath info as read-only
        with open(domTblPath, "r") as domTblData:
            # skip three header lines
            for i in range(3):
                domTblData.readline()

            dataList = domTblData.split()

            # since we split domTblData along whitespace, we need to reconstruct
            # the description string at the end of the table line. This string
            # always begins at index 27 and runs until the end of the list
            descString = " ".join(dataList[27:])

            infoDict["name"] = dataList[0]
            infoDict["tlen"] = dataList[2]
            infoDict['eValue'] = dataList[7]
            infoDict['desc'] = descString

        with open(hmmStatPath) as hmmStatData:
            # skip header lines
            for i in range(10):
                hmmStatData.readline()

            dataList = hmmStatData.split()


parser = argparse.ArgumentParser()
parser.add_argument("inputDir", help="Path to directory of nucleotide count files. Expected format is [prophageName]Chart.txt")
parser.add_argument("outputDir", help="path to directory to store output .png files")
parser.add_argument("force", help="If output files already exist, overwrite them", action="store_true")
args = parser.parse_args()

inputDir = args.inputDir

# credit to pycruft in https://stackoverflow.com/questions/3207219/how-do-i-list-all-files-of-a-directory for code for grabbing file paths
filePaths = []
for (dirpath, dirnames, files) in walk(inputDir):
    filePaths.extend(files)

filePaths.sort()

# credit to Tim Anderson for advice on how to open files:
# runs through all paths in the directory, opening as read only
for filePath in filePaths:
    nucleotideList = []

    with open("%s/%s" % (inputDir, filePath), "r") as fileData:
        # skip header line
        fileData.readline()
        for line in fileData:
            line.rstrip("\n")
            nucleotideList.append(int(line))

    regMatch = re.match(r'(.+?)Chart\.txt', filePath)
    # extract name from file path
    prophageName = regMatch.group(1)

    #print(len(nucleotideList))


    plotList(list=nucleotideList, prophageName=prophageName, outputDir=args.outputDir)
