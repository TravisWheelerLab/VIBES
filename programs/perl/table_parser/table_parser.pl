#!/usr/bin/env perl

use strict;
use warnings;
# import classes
use ViralSeq;

my $seq = ViralSeq->new(
    name        => "test",
    aliSt       => "5",
    aliEn       => "3",
    fullSeq     => "1",
    negStrand   => 0,
);

print "Length: " .$seq->printStats(). "\n";
