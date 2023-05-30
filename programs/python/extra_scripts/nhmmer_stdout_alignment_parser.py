# Credit to Daphne Demekas for the alignment parsing code!
# Big idea: Write a Python script that takes in a list of genomes, aligns every genome to every other genome with
# PyHMMER, grabs the stdout HMMER output alignment, and uses it to determine percent ID (or some other relatedness stat)
import argparse
import os
import sys


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
            ali_info_list = ''.join(target_split).split("\n")[3].split()
            for target_text in target_split:
                target = target_text.split()[0]
                if target == query:
                    continue
                alignment_text = target_text.split("Alignment:")[1:]
                assert len(alignment_text) == 1
                lines = [l for l in alignment_text[0].split("\n") if len(l) > 0]
                score = ali_info_list[1]
                query_from = ali_info_list[4]
                query_to = ali_info_list[5]
                target_from = ali_info_list[7]
                target_to = ali_info_list[8]
                ali_info = f">{query}[{query_from}:{query_to}]--{target}[{target_from}:{target_to}]\n"
                if save_scores_separate:
                    filename = f"{target}-{query}_{score}.txt"
                    alignmentfile = open(f"{data_dir}/{filename}", "w")
                    alignmentfile.write(ali_info)
                    for line in lines[1:]:
                        if query in line:
                            line_split = line.split()
                            querysub = ''.join([t for t in line_split[2:-1]])
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

                    alignmentfile.write(ali_info)

                    for line in lines[1:]:
                        if query in line:
                            line_split = line.split()
                            querysub = ''.join([t for t in line_split[2:-1]])
                            alignmentfile.write(querysub + "\n")
                        elif target in line:
                            line_split = line.split()
                            targetsub = ''.join([t for t in line_split[2:-1]])
                            alignmentfile.write(targetsub + "\n" + "\n")
                alignmentfile.close()


def parse_args(sys_args: list) -> argparse.Namespace:
    parser = argparse.ArgumentParser(sys_args, description="Parse nhmmer stdout to capture alignments between target "
                                                           "and query sequences. For each stdout file, produces a "
                                                           "directory named after the input text file containing "
                                                           "files with each alignment and the alignment score")
    parser.add_argument("stdout_dir", type=str, help="Path to directory populated with nhmmer stdout in text files")
    parser.add_argument("alignment_dir", type=str, help="Path to parent output directory, where alignments will be "
                                                        "stored in subdirectories")
    parser.add_argument("--suppress_scores", help="Do not include scores in alignent ", action="store_true",
                        default=False)

    return parser.parse_args()


def _main():
    # get arguments from argparse
    args = parse_args(sys.argv[1:])
    stdout_dir = args.stdout_dir
    alignment_dir = args.alignment_dir
    suppress_score = args.suppress_scores

    # flip suppress_score: if it's True, we want to input False and vice versa
    parse_nhmmer_stdout_alis(stdout_dir, alignment_dir, not suppress_score)


if __name__ == "__main__":
    _main()