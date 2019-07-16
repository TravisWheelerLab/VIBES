import subprocess
import argparse
from os import walk, path
import re
import sys


def DomTblGen(prophageName, inputDir, hmmPath, outputDir):

    outputFile = "%s/%s.domtbl" % (outputDir, prophageName)
    inputFile = "%s/%s.fasta" % (inputDir, prophageName)

    # if file exists already, and --force wasn't used, print a warning and exit
    if path.isfile(outputFile) and not force:
        print("\nWarning: File %s already exists. If you wish to overwrite files, use --force" % (outputFile))
        exit()

    else:
        # run hmmscant to create .domtbl of prophage
        command = "hmmscant --domtblout %s %s %s" % (outputFile, hmmPath, inputFile)
        subprocess.run(command.split())

# credit to Viktor Kerkez in: https://stackoverflow.com/questions/18160078/how-do-you-write-tests-for-the-argparse-portion-of-a-python-module
def parseArgs(sysArgs):
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta_dir", help="Path to directory of .fasta files. Expected format is [sequenceName].fasta")
    parser.add_argument("hmm_db", help="Path to HMM database to run against .fasta files")
    parser.add_argument("output_dir", help="Path to directory to store output .domtbl files")
    parser.add_argument("--force", help="If output files already exist, overwrite them", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    args = parseArgs(sys.argv[1:])
    inputDir = args.fasta_dir
    hmmPath = args.hmm_db
    outputDir = args.output_dir
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

        DomTblGen(prophageName, inputDir, hmmPath, outputDir)
