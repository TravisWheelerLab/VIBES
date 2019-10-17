#!/usr/bin/env perl

# goes through genomes in a directory, searching for tRNAs with tRNAscan-SE. Uses esl-sfetch to extract tRNA sequences from genomes and save them
use strict;
use warnings;

use Getopt::Long;

# Command-line options
my $inputDir = '';
my $inputFasta = '';
my $outputDir = '';
my $prefix = '';
my $help = '';
my $verbose = 0;
my $force = '';

GetOptions (
    "input_dir=s"   => \$inputDir,
    "input_fasta=s" => \$inputFasta,
    "output_dir=s"  => \$outputDir,
    "prefix=s"      => \$prefix,
    "help"          => \$help,
    "verbose"       => \$verbose,
    "force"         => \$force
    )
or die("Unknown option, try --help\n");

if ($help) {
    help();
}

# store all files in input directory in array. We expect that all files in input dir are genomes
my @genomes = glob "$inputDir/*";
# create hash table to store tRNAs and their respective sequences
my %trnaSeq;

if ($inputDir && $inputFasta) {
    die("--input_dir and --input_fasta can't be specificed at the same time!\n");
}

if ($inputDir) {
    # run tRNAscan-SE to search for tRNAs in genome
    foreach my $genome (@genomes) { 
        scanGenome($genome, %trnaSeq);
    }
}
else {
    scanGenome($inputFasta, %trnaSeq);
}

# if prefix has been specified, stick an '_' character on the end
if ($prefix) {
    $prefix = "$prefix\_";
}

foreach my $aminoAcid (keys %trnaSeq) {
    my $outputPath = "$outputDir/$prefix$aminoAcid.fasta";

    if (-e $outputPath) {
        unless($force) {
            die "$outputPath already exists. Use --force to automatically overwrite it";
        }
    }

    open(my $outputHandle, ">", $outputPath) or die "Can't open output file: $!\n";
    print $outputHandle $trnaSeq{$aminoAcid};
    close($outputHandle);
}

sub do_cmd {
    my $cmd = $_[0];

    if ($verbose > 0) {
        print "$cmd\n";
    }
    my $res = `$cmd 2>&1`;

    if ($verbose > 0) {
        print "$res\n";
    }
    
    return $res;
}

# scan genome for tRNAs and store them in trnaSeq
sub scanGenome {
    my $genome = shift;
    my %trnaSeq = shift;

    # run tRNAscan-SE with brief and quiet, which suppress header lines in output. Run with -B flag to specify bacterial tRNA
    my $tRNAscanResults = do_cmd("tRNAscan-SE -B --brief --quiet $genome");
    # index genome so esl-sfetch can search it
    do_cmd("esl-sfetch --index $genome");

    foreach my $line (split /\n/, $tRNAscanResults) {
        # use regex to capture info we want from tRNAscan-Se output line
        if ($line =~ /(.+?)\s+\d+\s+(\d+)\s+(\d+)\s+(\w+)\s+(.+?)\s+/) {
            my $seqName = $1;
            my $startCoord = $2;
            my $endCoord = $3;
            my $aminoAcid = $4;
            my $antiCodon = $5;

            # match as many characters as possible until we encounter the last '/', then grab everything until we come across a '.' character.
            # This should capture the name of the genome from its path
            $genome =~ /.+\/(.+?)\./;
            my $genomeName = $1;

            # esl-sfetch will output a .fasta format entry of the tRNA seq. We provide it with entry header
            my $fastaHeader = "$aminoAcid\_$genomeName\_$startCoord\.\.$endCoord Anticodon: $antiCodon";
            my $seqFasta = do_cmd("esl-sfetch -n \"$fastaHeader\" -c $startCoord..$endCoord $genome \"$seqName\"");

            # check if hash entry for this amino acid already exists. If so, add .fasta entry to it. Else, create entry and populate with .fasta info
            if (exists $trnaSeq{$aminoAcid}) {
                $trnaSeq{$aminoAcid} .= "$seqFasta\n";
            }
            else {
                $trnaSeq{$aminoAcid} = "$seqFasta\n";
            }

        }
        else {
            die "Regex failed to parse tRNAscan-SE output line $line\n";
        }
    }
    # clean up index file once we're done with it
    do_cmd("rm $genome.ssi");
}

sub help {
    print("
# find_trna.pl: Automatically runs tRNAscan-SE on all genomes in a directory. 
#The program uses the -B flag to search for bacterial tRNAs.
    Usage: perl find_trna.pl --input_dir /path/ --output_dir /path/

    Input:
        --input_dir /path/: Path to directory of target genomes in .fasta
            format. All files in this directory are expected to be genomes.
            This option is mutually exclusion with --input_fasta.

    Output:
        --output_dir /path/: Path to directory where output tRNA sequences
            will be stored in .fasta format. Since these files are intended
            to generate a MSA, tRNA sequences are saved in files named after
            the respective amino acids they encode for.
        --prefix str: Prefix that will be appended to output file names.

    Misc:
        --help: Prints this help page.
        --force: Overwrite output files, if they already exist.
        --verbose: Prints commands run by script and results of those");
}