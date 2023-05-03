# Credit to Daphne Demekas for the alignment parsing code!
# Big idea: Write a Python script that takes in a list of genomes, aligns every genome to every other genome with
# PyHMMER, grabs the stdout HMMER output alignment, and uses it to determine percent ID (or some other relatedness stat)
import subprocess
import os
import argparse


def parse_nhmmer_stdout_alis(directory, save_dir, save_scores_separate=True):
    """Given a directory containing text files with nhmmer output,
    this will find the alignments in the files
    and save them to the savedir

    The files will be saved as follows

    for each file in directory F, there will be a directory in savedir
    F_alignments

    Within each directory, there will be one file per target-query pair
    this """
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)
    for inputfile in os.listdir(directory):
        file = open(f'{directory}/{inputfile}', "r")
        alltext = file.read()

        query_split = alltext.split("Query:")[1:]

        data_dir = f"{save_dir}/{inputfile[:-4]}_alignments"
        if not os.path.exists(data_dir):
            os.mkdir(data_dir)
        for query_text in query_split:
            query = query_text.split()[0]
            target_split = query_text.split(">>")[1:]
            for target_text in target_split:
                target = target_text.split()[0]
                if target == query:
                    continue
                alignment_text = target_text.split("Alignment:")[1:]
                assert len(alignment_text) == 1
                lines = [l for l in alignment_text[0].split("\n") if len(l) > 0]
                lines_copy = lines.copy()
                if save_scores_separate:
                    score = lines_copy[0].strip().split()[1]
                    filename = f"{target}-{query}_{score}.txt"
                    alignmentfile = open(f"{data_dir}/{filename}", "w")
                    for line in lines[1:]:
                        if query in line:
                            line_split = line.split()
                            querysub = '>' + ''.join([t for t in line_split[2:-1]])
                            alignmentfile.write(querysub + "\n")
                        elif target in line:
                            line_split = line.split()
                            targetsub = ''.join([t for t in line_split[2:-1]])
                            alignmentfile.write(targetsub + "\n" + "\n")
                else:
                    filename = f"{target}-{query}.txt"
                    if os.path.exists(f"{data_dir}/{filename}"):
                        alignmentfile = open(f"{data_dir}/{filename}", "a")
                    else:
                        alignmentfile = open(f"{data_dir}/{filename}", "w")

                    for line in lines[1:]:
                        if query in line:
                            line_split = line.split()
                            querysub = '>' + ''.join([t for t in line_split[2:-1]])
                            alignmentfile.write(querysub + "\n")
                        elif target in line:
                            line_split = line.split()
                            targetsub = ''.join([t for t in line_split[2:-1]])
                            alignmentfile.write(targetsub + "\n" + "\n")
                alignmentfile.close()

def _main():
    test_stdout = "test/all_vs_all_alis/test_stdout.txt"
    test_alignment_dir = "test/all_vs_all_alis/alis/"

    parse_nhmmer_stdout_alis("test/all_vs_all_alis/stdout/", test_alignment_dir)


if __name__ == "__main__":
    _main()