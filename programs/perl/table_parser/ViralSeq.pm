#!/usr/bin/env perl
package ViralSeq;

{
use strict;
use warnings;
use diagnostics;

use Moose;

has name => (
    is          => 'ro',
    isa         => 'Str',
    required    => 1,
);

#First index of alignment
has aliSt => (
    is          => 'ro',
    isa         => 'Num',
    required    => 1,
);

#Last index of alignment
has aliEn => (
    is          => 'ro',
    isa         => 'Num',
    required    => 1,
);

#path to file containing complete viral sequence
has completeSeqPath => (
    is          => 'ro',
    isa         => 'Str',
    required    => 1,
);

#path to bacterial genome sequence was found in
has genomePath => (
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

has completeSeqLength => (
    is      => 'ro',
    lazy    => 1,
    builder => '_buildCompLength',
);

has percentComplete => (
    is      => 'ro',
    lazy    => 1,
    builder => '_buildPercent',
);

sub printStats {
    my $self = shift;
    my $returnString = '';

    return  "Length: " . $self->seqLength() . "\nComplete sequence length: "
            . $self->completeSeqLength() . "\nSize relative to complete sequence: "
            . $self->percentComplete();
}

sub _buildPercent {
    my $self = shift;
    return $self->seqLength() / $self->completeSeqLength();
}

sub _buildLength {
    my $self = shift;
    my $aliSt = $self->aliSt();
    my $aliEn = $self->aliEn();

    return abs($aliEn - $aliSt) + 1;
}

#Returns length of complete version of sequence
sub _buildCompLength {
    my $self = shift;
    open(my $fileHandle, "<", $self->completeSeqPath()) or die "Can't open .fna file " .$self->completeSeqPath() . ": $!";
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
