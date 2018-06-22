#!/usr/bin/env perl
package ViralSeq;

{
use strict;
use warnings;

use Moose;

my $isFullCutoff = 50;

has name => (
    is          => 'ro',
    isa         => 'Str',
    required    => 1,
);

#Start index in genome
has gnSt => (
    is          => 'ro',
    isa         => 'Num',
    required    => 1,
);

#End index in genome
has gnEn => (
    is          => 'ro',
    isa         => 'Num',
    required    => 1,
);

#Path to genome sequence was extracted from
has genomePath => (
    is  => 'ro',
    isa => 'Str',
);

#Start index of hit in reference sequence
has refSt => (
    is          => 'ro',
    isa         => 'Num',
    required    => 1,
);

#End index of hit in reference sequence
has refEn => (
    is          => 'ro',
    isa         => 'Num',
    required    => 1,
);

#path to file containing complete viral sequence
has referenceSeqPath => (
    is          => 'ro',
    isa         => 'Str',
    required    => 1,
);

#Is sequence on positive strand?
has isPos => (
    is          => 'ro',
    isa         => 'Bool',
    required    => 1,
);

has seqLength => (
    is      => 'ro',
    lazy    => 1,
    builder => '_buildLength',
);

has referenceSeqLength => (
    is      => 'ro',
    lazy    => 1,
    builder => '_buildRefLength',
);

has percentComplete => (
    is      => 'ro',
    lazy    => 1,
    builder => '_buildPercent',
);

has isFullLength => (
    is      => 'ro',
    lazy    => 1,
    builder => '_buildIsFull',
);

sub findFlankingAtts {
    my $self = shift;
    my $name = $self->name;
    my $genomePath = $self->genomePath;
    my $startIndex = $self->gnSt;
    my $endIndex = $self->gnEn;
    my $flankSize = 500;

    open(my $genome, "<", $genomePath) or die "Can't open $genomePath: $!";

    my $fastaHeader = readline($genome);

    my $identifier;

    if ($fastaHeader =~ />(.+?) /) {
        $identifier = $1;
    }
    else {
        die "\nUnable to parse .fna header line with regex!\n";
    }

#Five and three refer to the 5` and 3` ends of DNA, which runs from 5` to 3`
    my $fiveBegin;
    my $fiveEnd;
    my $threeBegin;
    my $threeEnd;

    if ($self->isPos) {
        #Since nucleotide indexing begins at 1, we want these at least == 1
        unless (($startIndex - $flankSize) < 1) {
            $fiveBegin = $startIndex - $flankSize;
        }
        else {
            $fiveBegin = 1;
        }

        unless (($startIndex - 1) < 1) {
            $fiveEnd = $startIndex - 1;
        }
        else {
            $fiveEnd = 1;
        }

        #These occuring beyond end of genome are handled later
        $threeBegin = $endIndex + 1;
        $threeEnd = $endIndex + $flankSize;
    }
    else { #is on negative strand
        #These occuring beyond end of genome are handled later
        $fiveBegin = $startIndex + $flankSize;
        $fiveEnd = $startIndex + 1;

        #Since nucleotide indexing begins at 1, we want these at least == 1
        unless (($endIndex - 1) < 1) {
            $threeBegin = $endIndex - 1;
        }
        else {
            $threeBegin = 1;
        }

        unless (($endIndex - $flankSize) < 1) {
            $threeEnd = $endIndex - $flankSize;
        }
        else {
            $threeEnd = 1;
        }
    }

#HEY ME: Split into 2 statements, check > 0
    do {

        my $fiveOutput = `esl-sfetch -c $fiveBegin..$fiveEnd $genomePath $identifier > fivePrimeFlank.fasta`;

    }while ();

    do {
        my $threeOutput = `esl-sfetch -c $threeBegin..$threeEnd $genomePath $identifier > threePrimeFlank.fasta`;
    }while ();

    `nhmmer -f $outputPath fivePrimeFlank.fasta threePrimeFlank.fasta`;

    `rm fivePrimeFlank.fasta threePrimeFlank.fasta`

}

sub tableLine {
    my $self = shift;

    my $returnString =  $self->name . "\t" . $self->gnSt . "\t" . $self->gnEn
                        . "\t" . $self->isFullLength . "\t" . $self->refSt
                        . "\t" . $self->refEn . "\t" . $self->isPos
                        . "\t" . $self->referenceSeqPath;

    if ($self->genomePath) {
        $returnString .= "\t" . $self->genomePath;
    }

    return  $returnString;
}

sub tableHeader {
    return "Name\tGnSt\tGnEn\tIsFullLength\tRefSt\tRefEn\tIsPos\tRefFile";
}

sub _buildIsFull {
    my $self = shift;
    my $beginMatches = 0;
    my $endMatches = 0;

    if ($self->refSt <= $isFullCutoff) {
        $beginMatches = 1;
    }
    if (($self->referenceSeqLength - $self->refEn) <= $isFullCutoff) {
        $endMatches = 1;
    }

    if ($beginMatches && $endMatches) {
        return 1;
    }
    else {
        return 0;
    }
}

sub _buildPercent {
    my $self = shift;
    return $self->seqLength / $self->referenceSeqLength;
}

sub _buildLength {
    my $self = shift;
    my $gnSt = $self->gnSt;
    my $gnEn = $self->gnEn;

    return abs($gnEn - $gnSt) + 1;

#Returns length of complete version of sequence
sub _buildRefLength {
    my $self = shift;
    open(my $fileHandle, "<", $self->referenceSeqPath) or die "Can't open .fna file " .$self->referenceSeqPath . ": $!";
    readline($fileHandle); #skip header line

    my $fastaBody = "";

    while (my $line = <$fileHandle>) {
        $fastaBody = $fastaBody . $line;
    }

    return length($fastaBody);
}
}
no Moose;
__PACKAGE__->meta->make_immutable;
1; #obligatory EOF true value
