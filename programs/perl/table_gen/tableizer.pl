#!/usr/bin/env perl
#Takes in a path to a genome directoy, a path to a file containing reference
#.hmm, paths to output directories, a suffix to append to the file name,
#and a boolean that sets the value of verbose.
#Returns a .dfam table for each genome in the directory.

#Credit to Kaitlin Carey for help in understanding how to run things on the cluster!

use strict;
use warnings;
use Time::HiRes;

#Genome number from server starts at one, so decrement to match array indexing
my $genomeNumber = $ARGV[0] - 1;
my $genomeDir = $ARGV[1];
my $referencePath = $ARGV[2];
my $tableDir = $ARGV[3];
my $scannedTableDir = $ARGV[4]
my $suffix = $ARGV[5];
my $verbose = $ARGV[6];

#we assume every entry in the directory is a genome we're interested in
my @genomes = glob "$genomeDir/*";

#Extract file name from path, then run nhmmscan and dfamscan
my $genome = $genomes[$genomeNumber];
my @splitLine = split('/', $genome);
my $fileName;

if ($splitLine[-1] =~ /(.+?)\./) {
    $fileName = $1;
}
else {
    be_verbose();
    die "Unable to extract file name from path: $genome\n";
}

my $tablePath = "$tableDir/$fileName" . "_$suffix.dfam";
my $scannedPath = "$scannedTableDir/$fileName" . "_$suffix" . "_scanned.dfam";

my $stTime = [Time::HiRes::gettimeofday()];
do_cmd("nhmmscan --cpu 1 --dfamtblout $tablePath $referencePath $genome");
my $elapsed = Time::HiRes::tv_interval($stTime);
do_cmd("perl dfamscan.pl --dfam_infile  $tablePath --dfam_outfile $scannedPath");

if ($verbose > 0) {
    be_verbose();
}

sub be_verbose {
    my $path = `pwd`;

    print "Number: $genomeNumber\n";
    print STDERR "Number: $genomeNumber\n";
    print "GenDir: $genomeDir\n";
    print "Ref: $referencePath\n";
    print "Table dir: $tableDir\n";
    print "Suffix: $suffix\n";
    print "Verbose: $verbose\n";
    print "Genome: $genome\n";
    print "Table Path: $tablePath\n";
    print "Scanned table path: $scannedPath\n";
    print "Time to complete nhmmscan: $elapsed seconds\n";
    print "Current directory: $path";
}

sub do_cmd {
    my $cmd = $_[0];

    if ($verbose > 0) {
        print "$cmd\n";
    }

    return `$cmd`;
}
