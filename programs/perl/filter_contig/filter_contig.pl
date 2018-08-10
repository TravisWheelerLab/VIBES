#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $input = '';
my $output = '';
my $cutoff = 1000;
my $help = 0;

my @contigArray = ();

GetOptions (
    "input=s"   => \$input,
    "output=s"  => \$output,
    "cutoff=i"  => \$cutoff,

    "help"      => \$help
    )
or die("Unknown option, try --help\n");

if ($help) {
    help();
    exit;
}

open(my $fileHandle, "<", $input) or die "Can't open input .fasta file: $!";
{
        local $/ = ">"; # use > instead of /n to end line
        readline $fileHandle; # skip first > character

        while (my $contig = <$fileHandle>) { #iterate through entries in .fasta

            #grab sequence, not header line, while leaving $contig intact
            if ($contig =~ /\d\n([ACTG\n]+)\n/) {
                my $sequence = $1;
                #remove any whitespace or > characters
                $sequence =~ s/\>|\s//g;
                my $length = length($sequence);

                if ($cutoff <= $length) {
                    #trim off >
                    $contig =~ s/\>//g;
                    push @contigArray, $contig;
                    }
            }
            else {
                die "Contig didn't match regex statement: $!";
            }

        }
}

open(my $outputHandle, ">", $output) or die "Can't open output .fasta file: $!";

foreach my $entry (@contigArray) {
    print $outputHandle ">$entry\n";
}

sub help {
    print "
    #contig_filter.pl: Takes in a .fasta file of contigs and filters them based
    #on length. All contigs at or above the cutoff length will be included in
    #the output file.

    Options:
        --input: Input .fasta file of unfiltered contigs
        --output: Output, filtered .fasta file. Will be overwritten by program
        --cutoff: Minumum length of contigs allowed by filter
        --help: Displays this help screen\n
    ";
}
