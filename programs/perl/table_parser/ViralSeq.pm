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

    open(my $FivePrime, ">", "./fivePrime.fasta") or die "Can't open fivePrime.fasta: $!";
    open(my $ThreePrime, ">", "./threePrime.fasta") or die "Can't open threePrime.fasta: $!";
    open(my $genome, "<", $genomePath) or die "Can't open $genomePath: $!";

    my $fastaHeader = readline($genome);

    my $identifier;

    if ($fastaHeader =~ />(.+?) /) {
        $identifier = $1;
    }
    else {
        die "\nUnable to parse .fna header line with regex!\n";
    }

    my $fiveBegin;
    my $fiveEnd;
    my $threeBegin;
    my $threeEnd;

    if ($self->isPos) {

    }
    else {

    }

    print $fivePrime `esl-sfetch -c $fiveBegin..$fiveEnd $genomePath $identifier`;
    print $threePrime `esl-sfetch -c $threeBegin..$threeEnd $genomePath $identifier`;

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
