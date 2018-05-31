#!/usr/bin/env perl

use strict;
use warnings;
# import classes
use ViralSeq;

#Path to folder containing reference prophage .fastas
my $prophagePath = $ARGV[0];
my %chartHash;

open(my $tableFile, "<", $ARGV[1]) or die "Can't open $ARGV[1]: $!";

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
        my $negStrand = 0;

        if ($strand eq '-'){
            my $negStrand = 1;
        }

        my $seq = ViralSeq->new(
            name                => $name,
            refSt               => $2,
            refEn               => $3,
            referenceSeqPath    => "$prophagePath$name.fasta",
            negStrand           => $negStrand,
            gnSt                => $5,
            gnEn                => $6,
        );

        unless (exists $chartHash{$name}) {
            my @array = 0x($seq->referenceSeqLength());
            my $arrayRef = \@array;

            $chartHash{$name} = $arrayRef;
        }

        #Each array index corresponds to a nucleotide in the reference strain.
        #Iterate over nucleotides present in genomic strain, incrementing
        #for each present in genome. To match array 0-indexing, nucleotide
        #indexes are offset by 1.
        for (my $i = $self->refSt() - 1; $i < $seq->refEn(); $i++) {
            my $arrayRef = $chartHash{$name};

            $arrayRef->[$i] += 1;
        }
    }
}
