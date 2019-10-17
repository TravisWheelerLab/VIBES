#!/usr/bin/env perl

# runs mafft on all .fasta files in a directory, outputing resultant aligned fastas in output dir 
use strict;
use warnings;

use Getopt::Long;

# Command-line options
my $inputDir = '';
my $outputDir = '';
my $prefix = '';
my $help = '';
my $verbose = 0;
my $force = '';

GetOptions (
    "input_dir=s"   => \$inputDir,
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
my @fastas = glob "$inputDir/*";

foreach my $fasta (@fastas) {
    # match as many characters as possible until we encounter the last '/', then grab everything until we come across a '.' character.
    # This should capture the name of the genome from its path
    $fasta =~ /.+\/(.+?)\./;
    my $fastaName = $1;

    do_cmd("mafft --auto $fasta > $outputDir/$fastaName.afa");
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

sub help {
    print();
}