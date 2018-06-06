#!/usr/bin/env perl
#Takes in a path to a genome directoy, a path to a file containing reference
#.hmm, and a path to an output directory. Returns a .dfam table for each
#genome in the directory.

use strict;
use warnings;

my $genomeDir = $ARGV[0];
my $referencePath = $ARGV[1];
my $tableDir = $ARGV[2];

#we assume every entry in the directory is a genome we're interested in
my @genomes = glob "$genomeDir/*";

#Extract file name from path, then run nhmmscan and dfamscan
foreach my $genome (@genomes) {
    my @splitLine = split('/', $genome);
    my $fileName;

    if ($splitLine[-1] =~ /(.+?)\./) {
        $fileName = $1;
    }

    my $tablePath = "$tableDir/$fileName" . "_viral.dfam";
    my $scannedPath = "$tableDir/$fileName" . "_viral_scanned.dfam";

    print "nhmmscan --dfamtblout $tablePath $referencePath $genome\n";
    print "perl dfamscan.pl --dfam_infile  $tablePath --dfam_outfile $scannedPath";
}
