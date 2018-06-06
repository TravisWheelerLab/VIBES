#!/usr/bin/env perl
#Takes in a path to a genome directoy, a path to a file containing reference
#.hmm files, and a path to an output directory. Returns a .dfam table for each
#genome in the directory.
#Arg 0: genome directory

use strict;
use warnings;

my $dir = $ARGV[0];
#we assume every entry in the directory is a genome we're interested in
my @genomes = glob "$dir/*";
