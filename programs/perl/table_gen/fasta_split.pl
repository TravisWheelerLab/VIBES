#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $inputFasta = '';
my $outputDir = '';
my $help = 0;

GetOptions(
    "input=s"   => \$inputFasta,
    "output=s"  => \$outputDir,
    "help"      => \$help
)
or die("Unknown argument, try --help\n");

if ($help) {
    help();
    exit;
}

open my $FH, "<$inputFasta" or die "Can't open $inputFasta: $!"; 
{
    local $/ = "\n>"; # set end of input record separator to \n>, so it splits on lines starting with >. Unfortunately,
    # this is necessary because very rarely, .fasta header lines will contain > characters somewhere in the middle

    # grab each entry, save it in a temporary file, run it through hmmbuild, and delete it
    while (my $entry = <$FH>)  {
        # remove any > characters at beginning or end of entry, leaving any in the middle of header lines intact
        $entry =~ s/>$//;
        $entry =~ s/^>//;
        # capture header, sequence data separately
        $entry =~ m/(.+?)\n(.+)/s;
        my $header = $1;
        my $seq = $2;
        # remove any whitespace in sequence data
        $seq =~ s/\s//g;

        # grab header line up to first whitespace character
        $header =~ m/(\S+)/;
        my $name = $1;

        my $outputPath = "$outputDir/$name.fasta";

        open my $outputHandle, ">$outputPath" or die "Can't open $outputPath: $!";
        {
            print $outputHandle ">$header\n$seq";
        }
        close $outputHandle;
    }
}

close $FH;
