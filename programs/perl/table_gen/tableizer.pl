#!/usr/bin/env perl
#Takes in a path to a genome directoy, a path to a file containing reference
#.hmm, a path to an output directory, and a suffix to append to the file name.
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

`nhmmscan --dfamtblout $tablePath $referencePath $genome`;
`perl dfamscan.pl --dfam_infile  $tablePath --dfam_outfile $scannedPath`;
