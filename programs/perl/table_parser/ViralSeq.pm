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

has evalue => (
    is          => 'ro',
    isa         => 'Num',
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
has referenceSeqLength => (
    is          => 'ro',
    isa         => 'Num',
    required    => 1,
);

#Is sequence on positive strand?
has strand => (
    is          => 'ro',
    isa         => 'Str',
    required    => 1,
);

has verbose => (
    is          => 'ro',
    isa         => 'Bool',
    required    => 1,
);

has flankingAttDir => (
    is      => 'ro',
    isa     => 'Str',
    default => "./attSites",
    lazy    => 1,
);

has isPos => (
    is      => 'ro',
    isa     => 'Bool',
    lazy    => 1,
    builder => '_buildIsPos'
);

has seqLength => (
    is      => 'ro',
    lazy    => 1,
    builder => '_buildLength',
);

has bacGenomeLength => (
    is      => 'ro',
    lazy    => 1,
    builder => '_buildBacGenomeLength');

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

has isFlanked => (
    is      => 'ro',
    writer  => '_set_isFlanked',
    default => 'Att Scan not run',
);

sub findFlankingAtts {
    my $self = shift;
    my $name = $self->name;
    my $genomePath = $self->genomePath;
    my $startIndex = $self->gnSt;
    my $endIndex = $self->gnEn;
    my $flankSize = 500;
    my $genomeSize = 0;
    #Five and three refer to the 5` and 3` ends of DNA, which runs from 5` to 3`
    my $fiveFile = $genomePath . "fivePrimeFlank.fasta";
    my $threeFile = $genomePath . "threePrimeFlank.fasta";
    my $fiveBegin;
    my $fiveEnd;
    my $threeBegin;
    my $threeEnd;
    my $fiveOutput = " ";
    my $threeOutput = " ";

    #Since we expect the genome to be in .fasta format, we open it to grab the
    #header line and then extract the string of characters before the first space,
    #which is used by esl-sfetch as an identifier. We also record the size of
    #the genome.
    open(my $genome, "<", $genomePath) or die "Can't open $genomePath: $!";

    my $fastaHeader = readline($genome);

    my $identifier;

    if ($fastaHeader =~ />(.+?) /) {
        $identifier = $1;
        #replace pesky pipes with escaped pipes
        $identifier =~ s/\|/\\|/g;
    }
    else {
        die "\nUnable to parse .fna header line ($fastaHeader) with regex!\n";
    }

    while (my $line = <$genome>) {
        chomp $line;
        $genomeSize += length $line;
    }

    if ($self->isPos) {
        #Since nucleotide indexing begins at 1, we want these at least == 1
        if (($startIndex - $flankSize) < 1) {
            $fiveBegin = 1;
        }
        else {
            $fiveBegin = $startIndex - $flankSize;
        }

        $fiveEnd = $startIndex;
        $threeBegin = $endIndex;

        #Don't let threeEnd exceed genome size
        if (($endIndex + $flankSize) > $genomeSize) {
            $threeEnd = $genomeSize;
        }
        else {
            $threeEnd = $endIndex + $flankSize;
        }
    }
    else { #is on negative strand
        if (($startIndex + $flankSize) > $genomeSize) {
            $fiveBegin = $genomeSize;
        }
        else {
            $fiveBegin = $startIndex + $flankSize;
        }
        $fiveEnd = $startIndex;
        $threeBegin = $endIndex;

        if (($endIndex - $flankSize) < 1) {
            $threeEnd = 1;
        }
        else {
            $threeEnd = $endIndex - $flankSize;
        }
    }

    #index genome so it's usable by esl-sfetch
    $self->_do_cmd("esl-sfetch --index $genomePath");

    $fiveOutput = $self->_do_cmd("esl-sfetch -n fivePrimeFlank -c $fiveBegin..$fiveEnd $genomePath $identifier");

    $threeOutput = $self->_do_cmd("esl-sfetch -n threePrimeFlank -c $threeBegin..$threeEnd $genomePath $identifier");

    open(my $fiveHandle, ">", $fiveFile) or die "Can't open $fiveFile: $!";
    open(my $threeHandle, ">", $threeFile) or die "Can't open $threeFile: $!";

    print $fiveHandle $fiveOutput;
    print $threeHandle $threeOutput;

    my $hmmerOutput = $self->_do_cmd("nhmmer --dna $fiveFile $threeFile");

    #The line [No hits detected that satisfy reporting thresholds] indicates that
    #no matching sequences were found, so don't create output file. In this case,
    #the viral seq is marked as flanked
    unless ($hmmerOutput =~ /\[No hits detected that satisfy reporting thresholds\]/) {
        my @splitLine = split('/', $genomePath);
        my $genomeName;

        if ($splitLine[-1] =~ /(.+?)\./) {
            $genomeName = $1;
        }

        my $outputPath = $self->flankingAttDir . "/$genomeName.$name.$startIndex.afa";

        open(my $outputHandle, ">", $outputPath) or die "Can't open $outputPath: $!";
        print $outputHandle $hmmerOutput;

        $self->_set_isFlanked(1);
    }
    else {
        $self->_set_isFlanked(0);
    }

    #clean up files we don't need
    $self->_do_cmd("rm $genomePath.ssi");
    $self->_do_cmd("rm $fiveFile $threeFile");
}

sub tableLine {
    my $self = shift;
    my $genomeName;

    # extract genome name from .fasta file path
    $self->genomePath =~ /([^\/]+)\./;
    $genomeName = $1;

    my $returnString =           $self->name
                        . "\t" . $self->evalue
                        . "\t" . $self->isFullLength
                        . "\t" . $self->refSt
                        . "\t" . $self->refEn
                        . "\t" . $self->referenceSeqLength
                        . "\t" . $self->isFlanked
                        . "\t" . $genomeName
                        . "\t" . $self->gnSt
                        . "\t" . $self->gnEn
                        . "\t" . $self->bacGenomeLength
                        . "\t" . $self->strand;

    return  $returnString;
}

sub tableHeader {
    return "Name\tE-Value\tIs Full Length Insertion\tMatch Start Position on Viral Genome\tMatch End Position on Viral Genome\tViral Genome Length\tFlanking Att Sites\tBacterial Genome Name\tMatch Start Position on Bacterial Genome\tMatch End Position on Bacterial Genome\tBacterial Genome Length\tStrand";
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

#Uses esl-seqstat to determine bacterial genome length, returns length
sub _buildBacGenomeLength {
    my $self = shift;
    my $bacGenomePath = $self->genomePath;
    my $bacGenomeLength;

    my $seqstatOutput = `esl-seqstat $bacGenomePath`;
    
    #capture total genome size from seqstat output with regex
    $seqstatOutput =~ /Total # residues:\s+?(\d+)\s/ or die "esl-seqstat regex failed for output $seqstatOutput: $!";
    $bacGenomeLength = $1;

    return $bacGenomeLength;
}

sub _buildIsPos {
    my $self = shift;
    my $isPos = 0;

    if ($self->strand == '+') {
        $isPos = 1;
    }

    return $isPos;
}

#Run command
sub _do_cmd {
    my $self = shift;
    my $cmd = $_[0];

    if ($self->verbose) {
        print "$cmd\n\n";
    }

    return `$cmd`;
}
}
no Moose;
__PACKAGE__->meta->make_immutable;
1; #obligatory EOF true value
