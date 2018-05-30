#!/usr/bin/env perl

use strict;
use warnings;
# import classes
use ViralSeq;

my $seq = ViralSeq->new(
    name            => "test",
    aliSt           => 5,
    aliEn           => 3,
    completeSeqPath => "../../../sequence_files/prophage/Pf_prophage_423117-5435766_strain_M1608.fasta",
    genomePath      => "test.fna",
    negStrand       => 1,
    isNeg           => 1,
);

print $seq->printStats(). "\n";
