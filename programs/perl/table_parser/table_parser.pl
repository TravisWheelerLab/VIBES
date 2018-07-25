#!/usr/bin/env perl
#Arguments:
    #0: path to folder containing reference prophage sequences
    #1: path to DFAM table
    #2: Path to genome
    #3: path to output 'bed-like' file
    #4: Path to output chart directory
    #5: Path to output flanking att side directory
    #6: Boolean isVerbose (prints commands to stdOut)


use strict;
use warnings;
use FindBin;
use lib $FindBin::Bin;
# import classes
use ViralSeq;

#Path to folder containing reference prophage .fastas
my $prophagePath = $ARGV[0];
my $genomePath = $ARGV[2];
my $attSitePath = $ARGV[5];
my $isVerbose = $ARGV[6];
my %chartHash;
my $chartPath;

open(my $tableFile, "<", $ARGV[1]) or die "Can't open $ARGV[1]: $!";
open(my $outputTSF, ">", $ARGV[3]) or die "Can't open $ARGV[3]: $!";


my $count = 1;

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
            name                => $name . "_$count",
            refSt               => $2,
            refEn               => $3,
            referenceSeqPath    => "$prophagePath$name.fasta",
            isPos               => $isPos,
            gnSt                => $5,
            gnEn                => $6,
            attSitePath         => $attSitePath,
            genomePath          => $genomePath,
            verbose             => $isVerbose,
        );

        $chartPath = "$ARGV[4]$name" . "Chart.txt";

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
            my $size = $seq->referenceSeqLength;
            if ($i < 0 || $i >= $size) {
                die "Nucleotide index is outside of bounds 1-$size for $name";
            }
            my $arrayRef = $chartHash{$name};
            $arrayRef->[$i] += 1;
        }

        if ($count == 1) {
            #Print header line
            print $outputTSF $seq->tableHeader . "\n";
        }

        print $outputTSF $seq->tableLine . "\n";
        #HMMER flanking sites
        $seq->findFlankingAtts;
        $count++;
    }
}

close $tableFile;

#since count is set to 1 initially for naming purposes, decrement by 1 for
#reporting accuracy
$count = $count - 1;

my @hashKeys = keys %chartHash;
print   "\nNumber of prophage sequences detected in $ARGV[1]: $count\nOf 50 " .
        "reference sequences, " . scalar @hashKeys . " had hits in $ARGV[1]\n";

#Print out index-based 'charts' in tab-delimited format
foreach my $hashKey (@hashKeys) {
    $chartPath = "$ARGV[4]$hashKey" . "Chart.txt";
    open(my $chartOutput, ">", $chartPath) or die "Can't open $chartPath: $!";

    my @chartArray = @{$chartHash{$hashKey}};
    my $chartLine = "$hashKey\n";

    foreach my $entry (@chartArray) {
        $chartLine .= "$entry\n";
    }

    print $chartOutput "$chartLine";
    close $chartOutput;
}
