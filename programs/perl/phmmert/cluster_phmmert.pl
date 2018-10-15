#!/usr/bin/env perl

use strict;
use warnings;
use Time::HiRes;
use Getopt::Long;

#command line argument variables that can be provided by user
my $genomeNumber;
my $genomeDir = '';
my $referenceHMMPath = '';
my $outputDir  = '';
my $verbose = 0;
my $help = 0;

#other variables
my $outputPath = '';
my $genome = '';
my $elapsed = '';

GetOptions (
    "referencehmm=s"    => \$referenceHMMPath,
    "outputdir=s"       => \$outputDir,
    "jobnumber=i"       => \$$genomeNumber,
    "genomes=s"         => \$genomeDir,
    "verbose"           => \$verbose,
    "help"              => \$help
    )
or die("Unknown argument, try --help\n");

sub run_phmmer {
    #put all files in directory into an array
    my @genomes = glob "$genomeDir/*";
    $genome = $genomes[$genomeNumber];

    #Extract file name from genome path
    my @splitLine = split('/', $genome);
    my $fileName;

    if ($splitLine[-1] =~ /(.+?)\./) {
        $fileName = $1;
    }
    else {
        be_verbose();
        die "Unable to extract file name from path: $genome\n";
    }

    $outputPath = "$outputDir/$fileName.tbl";

    do_cmd("phmmert --tblout $outputPath $referenceHMMPath $genome");
}

sub do_cmd {
    my $cmd = $_[0];

    if ($verbose > 0) {
        print "$cmd\n";
    }

    return `$cmd`;
}

sub help {

}

sub be_verbose {
    my $path = `pwd`;

    print "Number: $genomeNumber\n";
    print STDERR "Number: $genomeNumber\n";
    print "GenDir: $genomeDir\n";
    print "RefHMM: $referenceHMMPath\n";
    print "Output dir: $outputDir\n";
    print "Verbose: $verbose\n";
    print "Genome: $genome\n";
    print "Time to complete nhmmscan: $elapsed seconds\n";
    print "Current directory: $path";
}