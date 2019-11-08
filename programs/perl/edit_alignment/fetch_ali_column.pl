#!/usr/bin/env perl

# runs hmmbuild on all .afa files in a directory, outputing resultant hmms in output dir 
use strict;
use warnings;

use Getopt::Long;

# Command-line options
my $coords = '';
my $inputPath = '';
my $outputPath = '';
my $help = 0;
my $force = '';

GetOptions (
    "coords=s"      => \$coords,
    "input=s"       => \$inputPath,
    "output=s"      => \$outputPath,
    "help"          => \$help,
    "force"         => \$force
    )
or die("Unknown option, try --help\n");

if ($help) {
    help();
}
else {
    open my $output, ">$outputPath" or die "Can't open $outputPath: $!\n";
    open my $input, "<$inputPath" or die "Can't open $inputPath: $!\n";
    local $/ = "\n>"; # set end of input record separator to \n>, so it splits on lines starting with >. Unfortunately,
    # this is necessary because very rarely, .fasta header lines will contain > characters somewhere in the middle

    # extract start, end coords from argument. Populate initial coordinate values with something that will break later
    # if regex doesn't work properly
    my $startCoord = -1;
    my $endCoord = -1;
    if ($coords =~ /(\d+)\.\.(\d+)/) {
        $startCoord = $1;
        $endCoord = $2;
    }
    else {
        die "Unable to parse --coords argument: Expect formart is ####..####\n";
    }

    # grab each entry, save it in a temporary file, run it through hmmbuild, and delete it
    while (my $entry = <$input>)  {
        # remove any > characters at beginning or end of entry, leaving any in the middle of header lines intact
        $entry =~ s/>$//;
        $entry =~ s/^>//;
        # capture header, sequence data separately
        $entry =~ m/(.+?)\n(.+)/s;
        my $header = $1;
        my $seq = $2;

        # remove any whitespace in sequence data
        $seq =~ s/[\s]//g;

        # add 1 to length to substring inclusive (we don't want to slice off the character at $endCoord)
        my $subSeqLength = $endCoord - $startCoord + 1;

        # offset startCoord by 1 to account for 0-indexing
        my $subSeq = substr($seq, $startCoord - 1, $subSeqLength);

        print $output ">$header\n";
        print $output "$subSeq\n";
    }

    close($input);
    close($output);
}