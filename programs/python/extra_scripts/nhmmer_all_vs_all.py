# Credit to Daphne Demekas for the PyHMMER and alignment parsing code!
from Bio import SearchIO

def main(query_fasta, target_fasta, outputfile, stdout_file, hmmer_max = False):
    if hmmer_max:
        cmd = f'phmmer --cpu 16 --max --tblout {outputfile} \
            "{query_fasta}" "{target_fasta}"'
    else:
        cmd = f'phmmer --cpu 16 --tblout {outputfile} \
            "{query_fasta}" "{target_fasta}"'
    hmmer = subprocess.run(cmd,
        shell=True,
        capture_output=True,
        check=True,
    )

    with open(stdout_file,"w") as file:
        file.write(hmmer.stdout.decode("utf-8"))
        file.close()


def parse_stdout(alignment_file_dir: str, stdout_path: str, write_full_seq=True, querysequences: dict = None,
                 targetsequences: dict = None):
    """Parse the saved off stdout files from running hmmer

    alignment_file_dir: the directory of the txt file where to save the alignment.
        the alignments will be saved as query_id-target_id

    stdout_path: the path of the hmmer stdout txt file that was saved off from running phmmer

    write_full_seq: This is True if you want to save the original full length sequence along with the alignment

    querysequences, targetsequences: dictionaries of type {id: sequence}. This is needed if write_full_seq = True

    The output files saved in alignment_file_dir will be of structure:

    >query_id & target_id
    aligned_seq1
    aligned_seq2
    full_seq1
    full_seq2
    """

    if write_full_seq:
        assert querysequences is not None, "No querysequences dictionary is given. This is needed to write full sequences"
        assert targetsequences is not None, "No targetsequences dictionary is given. This is needed to write full sequences"

    result = SearchIO.parse(stdout_path, "hmmer3-text")

    for qresult in tqdm.tqdm(result):
        query_id = qresult.id
        for hit in qresult:
            target_id = hit.id

            if hit.evalue > 1:
                continue

            alignment_file_path = f"{alignment_file_dir}/{query_id}-{target_id}.txt"

            with open(alignment_file_path, "w") as alignment_file:
                for hsp in hit:
                    alignments = hsp.aln
                    seq1 = str(alignments[0].seq)
                    seq2 = str(alignments[1].seq)
                    alignment_file.write(">" + query_id + " & " + target_id + "\n")
                    alignment_file.write(seq1 + "\n")
                    alignment_file.write(seq2 + "\n")
                    if write_full_seq:
                        fullseq1 = querysequences[query_id]
                        fullseq2 = targetsequences[target_id]
                        alignment_file.write(fullseq1 + "\n")
                        alignment_file.write(fullseq2 + "\n")