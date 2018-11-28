import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys
from os import walk

# read in all files in total counts dir via walk
# for each file, populate an array in which each index corresponds to a line
# plot as line graph


def plotList(list, filename):
    # Credit to https://matplotlib.org/gallery/lines_bars_and_markers/simple_plot.html#sphx-glr-gallery-lines-bars-and-markers-simple-plot-py
    # for structure of code below
    fig, ax = plt.subplots()
    ax.plot(list)

    ax.set(xlabel='Nucleotide Index', ylabel='Number Found', title=filename)
    ax.grid()

    # fig.savefig("test.png")
    plt.show()


mypath = sys.argv[1]
print(mypath)

# credit to pycruft in https://stackoverflow.com/questions/3207219/how-do-i-list-all-files-of-a-directory for code for grabbing file paths
filePaths = []
for (dirpath, dirnames, files) in walk(mypath):
    filePaths.extend(files)

filePaths.sort()

# credit to Tim Anderson for advice on how to open files
for fileName in filePaths:
    nucleotideList = []

    with open("%s/%s" % (mypath, fileName), "r") as fileData:
        # skip header line
        fileData.readline()
        for line in fileData:
            line.rstrip("\n")
            nucleotideList.append(int(line))

    plotList(nucleotideList, fileName)
