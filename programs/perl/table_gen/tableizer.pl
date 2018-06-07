#!/usr/bin/env perl
#Takes in a path to a genome directoy, a path to a file containing reference
#.hmm, a path to an output directory, and a suffix to append to the file name.
#Returns a .dfam table for each genome in the directory.

use strict;
use warnings;

my $genomeDir = $ARGV[0];
my $referencePath = $ARGV[1];
my $tableDir = $ARGV[2];
my $suffix = $ARGV[3];

unless ($suffix) {
    die "Please enter a valid suffix argument to help distinguish which reference"
        . " strains are being used. Order of arguments is: Path to genome dir,"
        . " path to reference strain file, path to output dir, suffix argument\n";
}

#we assume every entry in the directory is a genome we're interested in
my @genomes = glob "$genomeDir/*";

#Extract file name from path, then run nhmmscan and dfamscan
foreach my $genome (@genomes) {
    my @splitLine = split('/', $genome);
    my $fileName;

    if ($splitLine[-1] =~ /(.+?)\./) {
        $fileName = $1;
    }

    my $tablePath = "$tableDir/$fileName" . "_$suffix.dfam";
    my $scannedPath = "$tableDir/$fileName" . "_$suffix" . "_scanned.dfam";

    print "\nBeginning nhmmscan on $genome\n";
    `nhmmscan --dfamtblout $tablePath $referencePath $genome`;
    print "\nBeginnning dfamscan.pl on $tablePath\n";
    `perl dfamscan.pl --dfam_infile  $tablePath --dfam_outfile $scannedPath`;
}
