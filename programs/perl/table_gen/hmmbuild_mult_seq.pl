#!/usr/bin/env perl

use strict;
use warnings;

my $pathToFolder = $ARGV[0];
$pathToFolder =~ s/\w+.fasta//;
# open .fasta file, demarcate lines with > rather than /n
open(my $fileHandle, "<", $ARGV[0]) or die "Can't open .fasta file $ARGV[0]: $!";
{
    local $/ = ">"; # use > instead of /n
    readline $fileHandle; # skip first > character
    my $inc = 1;

    # grab each entry, save it in a file, and run it through hmmerbuild
    while (my $seqEntry = <$fileHandle>) {
        $seqEntry =~ s/>|\r|//g; # removes trailing > character, carriage returns
        my $seqEntry = ">" . $seqEntry; # adds leading > character

        my $singleFile = $seqEntry;
        my $strainName = ""; #used in progress print statement
        if ($singleFile =~ />(.+)\n|\r/) { #grab header line to use as file name
            $singleFile = $1;
            $strainName = $1;
        }
        else {
            die "Sequence entry didn't match regex statement";
        }

        $singleFile =~ s/\d - \d|\d – \d/-/g; #get rid of spaces between long and short dashes
        $singleFile =~ s/ |\(|\)/_/g; #replace spaces and parens with _
        $singleFile = "$pathToFolder$singleFile.fasta"; #append file path and type

        open(my $fastaFile, ">", $singleFile) or die "Can't create output file at $singleFile: $!";
        print $fastaFile "$seqEntry";

        my $outputFile = $singleFile;
        $outputFile =~ s/\.fasta/.hmm/; #Create .hmm with same name as .fasta file

        `hmmbuild $outputFile $singleFile`;
        print "hmmbuild round $inc: $strainName completed\n";
        $inc++;
    }
}
