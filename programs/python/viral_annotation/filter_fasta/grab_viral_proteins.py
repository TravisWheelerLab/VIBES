import re
import argparse
import sys


# containing argparse in a function can make testing functions easier
# credit to Viktor Kerkez in: https://stackoverflow.com/questions/18160078/how-do-you-write-tests-for-the-argparse-portion-of-a-python-module
def parseArgs(sysArgs):
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Path to unfiltered .fasta file containing both viral and non-viral protein sequence")
    parser.add_argument("output", help="Path to output filtered .fasta file containing viral protein sequence")
    parser.add_argument("-f", "--force", help="If output files already exist, overwrite them", action="store_true")
    parser.add_argument("-v", "--verbose", help="Print number of matching entries found", action="store_true")

    return parser.parse_args()


# Grab all entries from input .fasta that have 'virus' or 'phage' in their headers and save them in output .fasta file. Discard any that don't.
# Print how many matches were found
def filterEntries(inputFasta, outputFasta, force, verbose):
    # Counter to keep track of, and ultimate report the number of, matching entries found
    savedEntries = 0

    # Since this regex statement will be used many times, we compile it ahead of time to save a little time. This regex statement
    # looks for a '>' at the start of a line, followed by at least 0 non-newline characters, followed by either 'virus', 'viral', or 'phage'.
    # This is to ensure that matches are in header lines, rather than sequence data, since the statement ignores case
    virusRegex = re.compile(r'>[^\n\r]*(virus|viral|phage)', re.IGNORECASE)

    # This variable stores either 'w' or 'x'. 'w' tells it to overwrite the output if it already exists (--force set),
    # 'x' instructs open() to only open the output file if it doesn't already exist (--force not set),
    # credit for conditional assignment syntax to https://stackoverflow.com/questions/394809/does-python-have-a-ternary-conditional-operator
    openPermission = 'w' if force else 'x'

    try:
        with open(outputFasta, openPermission) as outputFileHandle:
            with open(inputFasta, "r") as inputFileHandle:
                inputData = inputFileHandle.read()
                # use regex to split all entries in the file along the > character, which denotes the beginning of a .fasta entry
                entriesList = re.split(r'>', inputData)

                for entry in entriesList:
                    # re-attach '>' to front of entry's header line
                    entry = ">%s" % entry

                    # if the header line contains virus or phage (regardless of upper or lower case), we want to save it to output.
                    # Ignore matches in any line after the header
                    if re.search(virusRegex, entry):
                        outputFileHandle.write(entry)
                        savedEntries += 1

        if verbose:
            print("%d entries contained virus, viral, or phage in their headers" % savedEntries)

    # catch case where --force isn't set and file already exists:
    except FileExistsError:
        print("Output file %s already exists! To overwrite it and its contents, run with --force" % outputFasta)


# check that we're running the program, rather than just grabbing its functions (perhaps for testing in another program)
if __name__ == "__main__":
    args = parseArgs(sys.argv[1:])

    filterEntries(args.input, args.output, args.force, args.verbose)
