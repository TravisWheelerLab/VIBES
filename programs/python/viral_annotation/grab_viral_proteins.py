import re
import argparse
import sys


# containing argparse in a function can make testing functions easier
# credit to Viktor Kerkez in: https://stackoverflow.com/questions/18160078/how-do-you-write-tests-for-the-argparse-portion-of-a-python-module
def parseArgs(sysArgs):
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Path to unfiltered .fasta file containing both viral and non-viral proteins")
    parser.add_argument("output", help="Path to output filtered .fasta file containing viral protein sequence")
    parser.add_argument("--force", help="If output files already exist, overwrite them", action="store_true")

    return parser.parse_args()


# Grab all entries from input .fasta that have 'virus' or 'phage' in their headers and save them in output .fasta file. Discard any that don't.
# Print how many matches were found
def filterEntries(inputFasta, outputFasta):
    # Boolean that we'll use to check if we want an entry's sequence data
    matchingEntry = False
    savedEntries = 0

    with open(outputFasta, "w") as outputFileHandle:
        with open(inputFasta, "r") as inputFileHandle:
            for line in inputFileHandle:
                # .fasta format has a header line denoted by >, followed by a line or lines of whatever sequence. If the file starts with >, it's
                # a header line
                if re.match(r'>', line):
                    # if the header line contains virus or phage (regardless of upper or lower case), we want to save it to output. Set boolean to
                    # true. Otherwise, set boolean to false and ignore
                    if re.search(r'(virus|phage)', line, re.IGNORECASE):
                        matchingEntry = True
                        outputFileHandle.write(line)
                        savedEntries += 1

                    else:
                        matchingEntry = False

                # else, it's a sequence line. If the entry matches, we want to save it; otherwise, we ignore it
                else:
                    if matchingEntry:
                        outputFileHandle.write(line)

    print("%d entries contained virus or phage in their headers" % savedEntries)


# check that we're running the program, rather than just grabbing its functions (perhaps for testing in another program)
if __name__ == "__main__":
    args = parseArgs(sys.argv[1:])

    filterEntries(args.input, args.output)
