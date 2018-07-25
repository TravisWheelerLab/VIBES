#!/usr/bin/env perl
#Argument:
#   0: Number used to grab DFAM table from directory
#   1: Reference prophage seq directory
#   2: DFAM table directory
#   3: Genome directory
#   4: Output .tsv file directory
#   5: Nucleotide chart directory
#   6: Flanking att site directory
#   7: Booelean to determine if verbose

use strict;
use warnings;

#Table number from server starts at one, so decrement to match array indexing
my $tableNumber = $ARGV[0] - 1;
my $prophagePath = $ARGV[1];
my $tableDir = $ARGV[2];
my $genomeDir = $ARGV[3];
my $tsvDir = $ARGV[4];
my $chartDir = $ARGV[5];
my $flankingAttDir = $ARGV[6];
my $isVerbose = $ARGV[7];
my $genome;
my $genomePath;
my $tsvPath;

my @tables = glob "$tableDir/*";

if ($tables[$tableNumber] =~ /([^\/]+)_viral_scanned\.dfam/) {
    $genome = $1;
}
else {
    die "Can't extract genome name from table path: $tables[$tableNumber] $!";
}

$tsvPath = "$tsvDir/$genome.tsv";
$genomePath = "$genomeDir/$genome.fna";


my $parserOutput = do_cmd("perl table_parser.pl $prophagePath $tables[$tableNumber] $genomePath $tsvPath $chartDir $flankingAttDir $isVerbose");
print "$parserOutput\n";

if ($isVerbose) {
    be_verbose();
}

sub do_cmd {
    my $cmd = $_[0];

    if ($isVerbose > 0) {
        print "$cmd\n";
    }

    return `$cmd`;
}

sub be_verbose {
    my $path = `pwd`;

    print "\nTable: $tables[$tableNumber]\n";
    print "Reference prophage sequence directory: $prophagePath\n";
    print "Genome directory: $genomeDir\n";
    print "Genome path: $genomePath\n";
    print "Genome Name: $genome\n";
    print "TSV Path: $tsvPath\n";
    print "Chart directory: $chartDir\n";
    print "Flanking att site directory: $flankingAttDir\n";
    print "Current directory: $path\n";
}
