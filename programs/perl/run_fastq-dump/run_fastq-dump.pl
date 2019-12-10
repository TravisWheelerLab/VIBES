#!/usr/bin/env perl

# takes in a file of Sequence Read Archive (SRA) accession numbers and an output folder. Runs fastq-dump on each acc number, placing resultant
# .fastq file in output folder
use strict;
use warnings;

use Getopt::Long;

# Command-line options
my $inputPath = '';
my $outputDir = '';
my $help = 0;
my $verbose = 0;
my $force = 0;

GetOptions (
    "input=s"       => \$inputPath,
    "output=s"      => \$outputDir,
    "help"          => \$help,
    "force"         => \$force,
    "verbose"       => \$verbose
    )
or die("Unknown option, try --help\n");

if ($help) {
    help();
    exit;
}
else {
    open(my $accessionFile, "<$inputPath");

    while (my $accession = <$accessionFile>) {
        # strip whitespace characters such as \n from accession number. 'g' flag instructs regex to keep going after first match
        $accession =~ s/\s//g;

        # unless an entry already exists for acession or --force was enabled
        unless(-e "$outputDir/$accession.fastq" || $force) {
            do_cmd("fastq-dump -O $outputDir $accession");
        }
        else {
            print("File already found for $accession\n");
        }
    }
}

sub help {

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