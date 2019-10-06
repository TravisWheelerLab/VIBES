#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $genomeDir = '';
my $outputDir = '';
my $cpu = 0;
my $verbose = 0;
my $help = 0;

GetOptions (
    "gen_dir=s"         => \$genomeDir,
    "output_dir=s"      => \$outputDir,
    "cpu=i"             => \$cpu,
    "verbose"           => \$verbose,
    "help"              => \$help
    )
or die("Unknown argument, try --help\n");

my @genomePaths = glob( $genomeDir . '/*' );

foreach my $genomePath (@genomePaths) {
    my $genomeName = '';
    if ($genomePath =~ /.+\/(.+?)./) {
        $genomeName = $1;
    }
    do_cmd("nhmmer --noali --dfamtblout $outputDir/$genomeName.dfam $genomePath $genomePath")
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
    
    return ;
}