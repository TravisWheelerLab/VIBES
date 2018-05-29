#!/usr/bin/env perl
package ViralSeq;

{
use strict;
use warnings;

use Moose;

has name => (
    is          => 'rw',
    isa         => 'Str',
    required    => 1,
);

#First index of alignment
has aliSt => (
    is          => 'rw',
    isa         => 'Num',
    required    => 1,
);

#Last index of alignment
has aliEn => (
    is          => 'rw',
    isa         => 'Num',
    required    => 1,
);

#path to file containing complete viral sequence
has fullSeq => (
    is          => 'rw',
    isa         => 'Str',
    required    => 1,
);

#Is sequence on negative strand?
has negStrand => (
    is          => 'rw',
    isa         => 'Bool',
    required    => 1,
);

sub printStats {
    my $self = shift;
    my $aliSt = $self->{aliSt};
    my $aliEn = $self->{aliEn};

    return abs($aliEn - $aliSt) + 1;
}
}
1; #obligatory EOF true value
