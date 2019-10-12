#!/usr/bin/env perl

# goes through genomes in a directory, searching for tRNAs with tRNAscan-SE. Uses esl-sfetch to extract tRNA sequences from genomes and save them
use strict;
use warnings;

use Getopt::Long;

# Command-line options
my $inputDir = '';
my $outputDir = '';
my $help = '';
my $verbose = 0;
my $force = '';

GetOptions (
    "input_dir=s"   => \$inputDir,
    "output_dir=s"  => \$outputDir,
    "help"          => \$help,
    "verbose"       => \$verbose,
    "force"         => \$force
    )
or die("Unknown option, try --help");

if ($help) {
    help();
}

my @genomes = glob "$inputDir/*";

# run tRNAscan-SE to search for tRNAs in genome
foreach my $genome (@genomes) { 
    # run tRNAscan-SE with brief and quiet, which suppress header lines in output
    my $tRNAscanResults = do_cmd("tRNAscan-SE --brief --quiet $genome");    
    foreach my $line (split /\n/, $tRNAscanResults) {
        if ($line =~ /(.+?)\s+\d+\s+(\d+)\s+(\d+)\s+(\w+)\s+([ACTG]+)/) {
            my $seqName = $1;
            my $startCoord = $2;
            my $endCoord = $3;
            my $aminoAcid = $4;
            my $antiCodon = $5;

            my $fastaHeader = 

            do_cmd("esl-sfetch --index $genome");
            my $seq = do_cmd("esl-sfetch -c $startCoord..$endCoord $genome \"$seqName\"");

        }
        else {
            die "Regex failed to parse tRNAscan-SE output line $line\n";
        }
    } 
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