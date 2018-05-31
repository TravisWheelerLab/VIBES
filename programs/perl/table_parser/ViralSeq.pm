#!/usr/bin/env perl
package ViralSeq;

{
use strict;
use warnings;
use diagnostics;

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

#Is sequence on negative strand?
has negStrand => (
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
    buidler => '_buildIsFull',
)

sub printTableLine {
    my $self = shift;
    my $returnString = '';

    return  "Length: " . $self->seqLength() . "\nReference sequence length: "
            . $self->referenceSeqLength() . "\nSize relative to reference sequence: "
            . $self->percentComplete();
}

sub _buildIsFull {
    my $self = shift;
    my $beginMatches = 0;
    my $endMatches = 0;

    if ($self->refSt() < $isFullCutoff) {
        $beginMatches = 1;
    }
    if (($self->referenceSeqLength() - $self->refEn()) < $isFullCutoff) {
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
    return $self->seqLength() / $self->referenceSeqLength();
}

sub _buildLength {
    my $self = shift;
    my $gnSt = $self->gnSt();
    my $gnEn = $self->gnEn();

    return abs($gnEn - $gnSt) + 1;
}

#Returns length of complete version of sequence
sub _buildRefLength {
    my $self = shift;
    open(my $fileHandle, "<", $self->referenceSeqPath()) or die "Can't open .fna file " .$self->referenceSeqPath() . ": $!";
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
