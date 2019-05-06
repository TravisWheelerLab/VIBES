import re
import argparse
import sys


# containing argparse in a function can make testing functions easier
# credit to Viktor Kerkez in: https://stackoverflow.com/questions/18160078/how-do-you-write-tests-for-the-argparse-portion-of-a-python-module
def parseArgs(sysArgs):
    parser = argparse.ArgumentParser()
    parser.add_argument("input.fasta", help="Path to unfiltered .fasta file containing both viral and non-viral proteins")
    parser.add_argument("output.fasta", help="Path to output filtered .fasta file containing viral protein sequence")
    parser.add_argument("--force", help="If output files already exist, overwrite them", action="store_true")


# Grab all entries from input .fasta that have 'virus' or 'phage' in their headers. Discard any that don't.
def filterEntries(inputFasta):

    matchResult = re.match(r'(virus|phage)', headerLine, re.IGNORECASE)


if __name__ == "__main__":
    print("Hi")
    args = parseArgs(sys.argv[1:])
