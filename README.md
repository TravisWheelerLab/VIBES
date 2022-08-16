# Pseudomonas Pipeline

This is a pipeline that annotates viral insertions in bacterial genomes, given a set of
reference viruses and bacterial genomes. Useful to bioinformaticians and computational biologists
because it can allow you to search for a large number of viruses in a large amount of bacterial
genomes with relatively little effort. The pipeline can be run in two ways: massively parallel
on a server cluster with a job manager like SLURM, or for small data sets, on a desktop or
laptop with GNU parallel.

## Pipeline

TODO: What does the pipeline look like? Can reference scripts described later,
bonus points if you include a diagram
The pipeline is a collection of scripts written in Python and Perl. Intended to be run as follows:

## Using

TODO: Answer the following questions in this section

  1. What data are required to run it and what format should they be in?
  2. Which versions of Perl and Python are supported?
  3. What dependencies are required for each script? Maybe add alongside
     descriptions
  4. What format will the ultimate output be in?

  ![Pipeline Diagram](pipeline_diagram-1-3.png)

  Input Data: 	- single .fasta file of reference viral genomes
  				- directory of reference bacterial genomes, one per .fasta file
  				- if annotation is desired: file of reference viral proteins, .hmm or .fasta format
  				- if annotation is desired: file of reference viral domains, .hmm or .fasta format

  Python Supported:  
  				- Python 3

  Output:  
  				- .json file for each bacterial genome describing detected
  				viral insertions  
				- optionally, .tsv (tab-separated value) file for each bacterial genome describing detected
				viral insertions  
  				- SODA plots showing location of each viral integration in bacterial genome

Running VIBES:  
				- Likely the easiest way to run VIBES is by running the Nextflow workflow, which accepts pipeline input,
				runs each step of the pipeline, automatically submits jobs to job managers such as SLURM, supports both
				Docker and Singularity, and will deposit output files in whatever directory is set in its .config file.
				See the README file in the nextflow_workflow/ dir for more details.  
				- If desired, users can also run each stage of the pipeline manually and submit manually-written
				massively parallel scripts to job managers. This is not the recommended method of running VIBES.

### Docker Image

There is a [Docker image](https://hub.docker.com/repository/docker/traviswheelerlab/pseudomonas_pipeline_runner)
capable of running the software, including all dependencies. As a convenience,
it can also be built using the `build-container.sh` script in the project root.

## Programs and Scripts

### hmmbuild_mult_seq.py
Dependencies: argparse, hmmbuild, hmmpress, os, re, subprocess, sys, typing

This program accepts an input .fasta file and converts its entries into HMMs using hmmbuild, storing them in the 
user-provided output 'pressed' .hmm file. User can specify program verbosity and how many threads hmmbuild will use. 
Users can specify input sequence type: dna, rna, or amino acid.

### tableizer.py
Dependencies: argparse, dfamscan.pl, mv, nhmmscan, os, subprocess, sys, typing

This program searches for viral genomes or genome fragments in input bacterial genomes using nhmmscan, saving results in
DFAM table format (nhmmscan's --dfamtblout option). Next, it automatically runs dfamscan.pl, a program that resolves
multiple competing matches to the same region of bacterial DNA by keeping only the best match. Required positional 
arguments are 1) hmm_db (.hmm file containing all viral genomes we want to search for), 2) genome_path (path to input 
bacterial genome), and 3) output_table (path to output scanned .dfam file). User can also specify program verbosity with
--verbose, the number of worker threads with --cpu (default 1), an optional path to dfamscan.pl with --dfamscan_path 
(default expects it to live in the same directory as this file), and whether to overwrite output files if they already 
exist with --force.

### table_parser.py
Dependencies: argparse, esl-seqstat, json, os, subprocess, sys, typing

This program parses scanned .dfam tables, filtering out integrations that fail to pass an e-value threshold. It
produces .tsv files describing detected integrations in a bacterial genome and .json files that contain reference viral
genome annotation information (if enabled) and an array in which each array index corresponds to a nucleotide index in
the reference viral genome the .json file is named after. The value of an array index represents how many integrations 
across all bacterial genomes in the present run of the pipeline included the nucleotide at the corresponding position in
the file's reference viral genome. Requires an input .dfam file (dfam_path), an input bacterial genome (genome_path), an
output .tsv file path (output_tsv_path), an output .json path (output_json_path). Optionally, the user can specify the
maximum allowable e-value for a potential integration with --max_evalue (default 1e-5, where lower is better), program
verbosity with --verbose, and whether output files that already exist should be overwritten with --force.

### grab_viral_proteins.py
Dependencies: re, argparse, sys

This program is only sort of part of the pipeline: it's designed to scan through
the headers of Uniprot protein .fasta files, looking for entries with keywords
"virus", "viral", or "phage". Produces an output .fasta file containing matching
fasta entries (both headers and sequence). This program was used to generate a
Uniprot viral protein database for viral annotation.

### generate_domtbls.py
Dependencies: hmmscant, subprocess, argparse, os, re, sys

This program uses hmmscant to search for viral proteins/domains contained in a
HMM database against .fasta format viral genomes, outputting results in .domtbl
format. Accepts directory of .fasta reference viral genomes, a file containing
viral sequence HMMs, an output dir, and a NCBI genetic code #.

### domtblscan.pl
Dependencies: nhmmscan, trf, Getopt::Long, File::Temp, Cwd, File::Basename, strict, warnings

Modified dfamscan.pl, changed to work on .domtbl format files. Adjudicates
between competing annotations at the same location in a genome

### annotation_methods.py
Dependencies: re

Describes Match and Genome objects, as well as a number of methods used to
annotate viral genomes with proteins and domains. The most important of these
methods is annotateGenomes, which returns a Genome object containing a
dictionary populated with Matches representing proteins/domains that matched to
the viral genome.

### visualize_nucleotide_counts.py
Dependencies: matplotlib, numpy, argparse, re, os.walk, sys, json, annotation_methods.py

Uses nucleotide index chart files and the annotateGenomes() method from
annotation_methods.py to draw a plot displaying the number of detected
occurences of each nucleotide in a viral genome. Below the x-axis,
differently-colored lines that correspond to a matching protein/domain are drawn
and labeled with their accession IDs. Because protein/domain matches can
overlap, lines are vertically staggered to prevent overlap. Accepts a directory
of nucleotide count charts, a directory of viral protein match .domtbl files, a
directory of domain match .domtbl files, a directory of .dfam files describing
matches to viral genes of particular interest to the virologists, and an ouput
directory where output .png plots are saved. Plots are drawn with matplotlib.

### integrase_dist_plot.py
Dependencies: re, argparse, os.walk, matplotlib, annotation_methods.py

Integrases seemed more likely to end up in bacterial genomes, so this program
goes counts how many integrase matches were detected per bacterial genome. It
does this by combining information about which regions of viruses are found in
bacterial genomes stored in bacterial .dfam files with information about what
proteins/domains are matched to regions of viral genomes found in viral
annotation .domtbl files. Any time a viral region corresponding to an integrase
is detected in a bacterial genome, its integrase count is iterated. This program
produces two types of output: a .csv file decribing how many integrases were
detected for every bacterial genome, and an output .pdf of a histrogram
describing the distribution of detected integrases.

## License

TODO(george): BSD 3-clause

## Authors

  - Conner Copeland
  - ...
  - Travis Wheeler
  - George Lesica

## Acknowledgements

???

