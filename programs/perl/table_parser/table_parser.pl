#!/usr/bin/env perl
#Arguments:
    #0: path to folder containing reference prophage sequences
    #1: path to DFAM table
    #2: path to output 'bed-like' file
    #3: Path to chart directory

use strict;
use warnings;
use FindBin;
use lib $FindBin::Bin;
# import classes
use ViralSeq;

#Path to folder containing reference prophage .fastas
my $prophagePath = $ARGV[0];
my %chartHash;

open(my $tableFile, "<", $ARGV[1]) or die "Can't open $ARGV[1]: $!";
open(my $outputTSF, ">", $ARGV[2]) or die "Can't open $ARGV[2]: $!";
my $chartPath;

my $count = 0;

while (my $line = <$tableFile>) {

    # Skip lines starting with #
    if ($line =~ /#/) {
    }

    #Well folks, it's already time for a monsterous, unreadable Regex line.
    #To summarize its twisted inner workings: Group 1 captures the name of the
    #strain, group 2 nabs the start index of the reference strain hit,
    #group 3 grabs the end index of the reference strain hit, group 4 obtains
    #the + or - sign indicating which strand the sequence is one, group 5 nicks
    #the start index in the bacterial genome, and group 6 contains the end index
    #in the bacterial genome.
    elsif ($line =~ /(.+?)\s+.+?\s+.+?\s+.+?\s+.+?\s+.+?\s+(.+?)\s+(.+?)\s+(.+?)\s+(.+?)\s+(.+?)\s/){
        my $name = $1;
        my $strand = $4;
        my $isPos = 0;

        if ($strand eq "+"){
            $isPos = 1;
        }

        my $seq = ViralSeq->new(
            name                => $name,
            refSt               => $2,
            refEn               => $3,
            referenceSeqPath    => "$prophagePath$name.fasta",
            isPos               => $isPos,
            gnSt                => $5,
            gnEn                => $6,
        );

        $chartPath = "$ARGV[3]$name.txt";

        unless (exists $chartHash{$name}) {
            my @array = (0)x$seq->referenceSeqLength;
            my $arrayRef = \@array;

            #if file already exists, read in its contents and populate array
            if (open(my $chartFile, "<", $chartPath)) {
                # skip header
                readline $chartFile;

                my $index = 0;
                while (my $line = <$chartFile>) {
                    chomp $line;
                    $arrayRef->[$index] = $line;
                    $index += 1;
                }

                close $chartFile;
            }
            $chartHash{$name} = $arrayRef;
        }

        #Each array index corresponds to a nucleotide in the reference strain.
        #Iterate over nucleotides present in genomic strain, incrementing
        #for each present in genome. To match array 0-indexing, nucleotide
        #indexes are offset by 1.
        for (my $i = $seq->refSt - 1; $i < $seq->refEn; $i++) {
            my $arrayRef = $chartHash{$name};
            $arrayRef->[$i] += 1;
        }

        if ($count == 0) {
            #Print header line
            print $outputTSF $seq->tableHeader . "\n";
        }

        print $outputTSF $seq->tableLine . "\n";
        $count++;
    }
}

my @hashKeys = keys %chartHash;
print   "\nNumber of prophage sequences detected in $ARGV[1]: $count\nOf 50 " .
        "reference sequences, " . scalar @hashKeys . " had hits in $ARGV[1]\n";

#Print out index-based 'charts' in tab-delimited format
foreach my $hashKey (@hashKeys) {
    $chartPath = "$ARGV[3]$hashKey" . "Chart.txt";
    open(my $chartOutput, ">", $chartPath) or die "Can't open $chartPath: $!";

    my @chartArray = @{$chartHash{$hashKey}};
    my $chartLine = "$hashKey\n";

    foreach my $entry (@chartArray) {
        $chartLine .= "$entry\n";
    }

    print $chartOutput "$chartLine";
    close $chartOutput
}
