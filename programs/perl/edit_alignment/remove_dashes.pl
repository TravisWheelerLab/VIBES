#!/usr/bin/env perl

# runs hmmbuild on all .afa files in a directory, outputing resultant hmms in output dir 
use strict;
use warnings;

use Getopt::Long;

# Command-line options
my $inputPath = '';
my $outputPath = '';
my $help = 0;
my $force = '';

GetOptions (
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
    open(my $output, ">", $outputPath) or die "Can't open $outputPath: $!\n";
    open(my $input, "<", $inputPath) or die "Can't open $inputPath: $!\n";

    while (my $line = <$input>) {
        # if line starts with >, we know it's a header line and we don't need to remove dash characters. Simply print to output
        # only check against ASCII characters, to avoid weird invisible unicode nonsense
        if ($line =~ /^>/a) {
            print $output "\n$line";
        }
        else {
            # use regex to remove all dash characters from body line
            $line =~ s/\-|\n//g;
            print $output $line;
        }
    }

    close($input);
    close($output);
}