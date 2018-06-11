#!/usr/bin/env perl
#Takes in a path to a genome directoy, a path to a file containing reference
#.hmm, a path to an output directory, a suffix to append to the file name,
#and a boolean that sets the value of verbose.
#Returns a .dfam table for each genome in the directory.

#Credit to Kaitlin Carey for help in understanding how to run things on the cluster!

use strict;
use warnings;

#Genome number is
my $genomeNumber = $ARGV[0];
my $genomeDir = $ARGV[1];
my $referencePath = $ARGV[2];
my $tableDir = $ARGV[3];
my $suffix = $ARGV[4];
my $verbose = $ARGV[5];

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
    die "Unable to extract file name from path: $genome\n";
}

my $tablePath = "$tableDir/$fileName" . "_$suffix.dfam";
my $scannedPath = "$tableDir/$fileName" . "_$suffix" . "_scanned.dfam";

do_cmd("nhmmscan --dfamtblout $tablePath $referencePath $genome");
do_cmd("perl dfamscan.pl --dfam_infile  $tablePath --dfam_outfile $scannedPath");

sub be_verbose {
    our $genomeNumber;
    our $genomeDir;
    our $referencePath;
    our $tableDir;
    our $suffix;
    our $genome;
    our $tablePath;
    our $scannedPath;
    my $path = `pwd`;

    print "Number: $genomeNumber\n";
    print STDERR "Number: $genomeNumber\n";
    print "GenDir: $genomeDir\n";
    print "Ref: $referencePath\n";
    print "Table dir: $tableDir\n";
    print "Suffix: $suffix\n";
    print "Genome: $genome\n";
    print "Table Path: $tablePath\n";
    print "Scanned table path: $scannedPath\n";
    print "Current directory: $path";
}

sub do_cmd {
    our $verbose;
    my $cmd = $_[0];

    if ($verbose > 0) {
        print "$cmd\n";
    }

    return `$cmd`;
}
