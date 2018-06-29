#!/usr/bin/env perl
package ViralSeq;

{
use strict;
use warnings;

use Moose;

my $isFullCutoff = 500;

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

#Path to genome the sequence was extracted from
has genomePath => (
    is          => 'ro',
    isa         => 'Str',
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

#Is sequence on positive strand?
has isPos => (
    is          => 'ro',
    isa         => 'Bool',
    required    => 1,
);

has verbose => (
    is          => 'ro',
    isa         => 'Bool',
    required    => 1,
);

has attSitePath => (
    is      => 'ro',
    isa     => 'Str',
    default => "./attSites",
    lazy    => 1,
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
    my $outputPath = $self->attSitePath . "/$name";
    my $flankSize = 2000;
    #Five and three refer to the 5` and 3` ends of DNA, which runs from 5` to 3`
    my $fiveBegin;
    my $fiveEnd;
    my $threeBegin;
    my $threeEnd;
    my $fiveOutput = " ";
    my $threeOutput = " ";

    open(my $genome, "<", $genomePath) or die "Can't open $genomePath: $!";

    my $fastaHeader = readline($genome);

    my $identifier;

    if ($fastaHeader =~ />(.+?) /) {
        $identifier = $1;
        #replace pesky pipes with escaped pipes
        $identifier =~ s/\|/\\|/g;
    }
    else {
        die "\nUnable to parse .fna header line with regex!\n";
    }

    if ($self->isPos) {
        #Since nucleotide indexing begins at 1, we want these at least == 1
        unless (($startIndex - $flankSize) < 1) {
            $fiveBegin = $startIndex - $flankSize;
        }
        else {
            $fiveBegin = 1;
        }

        $fiveEnd = $startIndex;

        #These occuring beyond end of genome are handled later
        $threeBegin = $endIndex;
        $threeEnd = $endIndex + $flankSize;
    }
    else { #is on negative strand
        #These occuring beyond end of genome are handled later
        $fiveBegin = $startIndex + $flankSize;
        $fiveEnd = $startIndex;

        $threeBegin = $endIndex;

        unless (($endIndex - $flankSize) < 1) {
            $threeEnd = $endIndex - $flankSize;
        }
        else {
            $threeEnd = 1;
        }
    }

    #index genome so it's usable by esl-sfetch
    $self->_do_cmd("esl-sfetch --index $genomePath");

    #Check that indexes do not exceed genome size using esl-sfetch
    do {
        if ($fiveOutput =~ /Subsequence end \d+ is greater than length (\d+)/){
            $fiveBegin = $1;
        }

        $fiveOutput = $self->_do_cmd("esl-sfetch -n $name -c $fiveBegin..$fiveEnd $genomePath $identifier");
    }while ($fiveOutput =~ /Subsequence end \d+ is greater than length (\d+)/);

    do {
        if ($threeOutput =~ /Subsequence end \d+ is greater than length (\d+)/){
            $threeEnd = $1;
        }

        $threeOutput = $self->_do_cmd("esl-sfetch -n $name -c $threeBegin..$threeEnd $genomePath $identifier");
    }while ($threeOutput =~ /Subsequence end \d+ is greater than length (\d+)/);

    open(my $fiveHandle, ">", "./fivePrimeFlank.fasta") or die "Can't open fivePrimeFlank.fasta: $!";
    open(my $threeHandle, ">", "./threePrimeFlank.fasta") or die "Can't open threePrimeFlank.fasta: $!";

    print $fiveHandle $fiveOutput;
    print $threeHandle $threeOutput;

    $self->_do_cmd("nhmmer -o $outputPath fivePrimeFlank.fasta threePrimeFlank.fasta");

    #clean up files we don't need
    $self->_do_cmd("rm $genomePath.ssi");
    #$self->_do_cmd("rm fivePrimeFlank.fasta threePrimeFlank.fasta");
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
}
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
#Run command
sub _do_cmd {
    my $self = shift;
    my $cmd = $_[0];

    if ($self->verbose > 0) {
        print "$cmd\n";
    }

    return `$cmd`;
}
}
no Moose;
__PACKAGE__->meta->make_immutable;
1; #obligatory EOF true value
