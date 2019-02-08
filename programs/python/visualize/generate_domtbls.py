import subprocess
import argparse
from os import walk, path
import re


def DomTblGen(prophageName, inputDir, outputDir):

    outputFile = "%s/%s.domtbl" % (outputDir, prophageName)
    inputFile = "%s/%s.fasta" % (inputDir, prophageName)

    # if file exists already, and --force wasn't used, print a warning and exit
    if path.isfile(outputFile) and not force:
        print("\nWarning: File %s already exists. If you wish to overwrite files, use --force" % (outputFile))
        exit()

    else:
        # run hmmscant to create .domtbl of prophage
        command = "hmmscant --domtblout %s ../../../sequence_files/hmm/TIGRFAM.HMM %s" % (outputFile, inputFile)
        subprocess.run(command.split())


parser = argparse.ArgumentParser()
parser.add_argument("inputDir", help="Path to directory of .fasta files. Expected format is [sequenceName].fasta")
parser.add_argument("outputDir", help="Path to directory to store output .domtbl files")
parser.add_argument("--force", help="If output files already exist, overwrite them", action="store_true")
args = parser.parse_args()

inputDir = args.inputDir
outputDir = args.outputDir
force = args.force

# credit to pycruft in https://stackoverflow.com/questions/3207219/how-do-i-list-all-files-of-a-directory for code for grabbing file paths
filePaths = []
for (dirpath, dirnames, files) in walk(inputDir):
    filePaths.extend(files)

filePaths.sort()

for filePath in filePaths:
    regMatch = re.match(r'(.+?)\.fasta', filePath)
    # extract name from file path
    prophageName = regMatch.group(1)

    DomTblGen(prophageName, inputDir, outputDir)
