#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib $FindBin::Bin;
use ViralSeq;
use Getopt::Long;

my $tablePath = "";
my $genomePath = "";
my $tsvDir = "";
my $chartDir = "";
my $flankingAttDir = "";
my $suffix = '';
my $verbose = ''; #default false value
my $help = ''; #^^
my $force = 0; #^^
my $dfamName;
my $tsvPath;
my $maxEval = 1e-5;

GetOptions (
    "dfam=s"        => \$tablePath,
    "genome=s"      => \$genomePath,
    "tsv=s"         => \$tsvDir,
    "index_charts=s"=> \$chartDir,
    "att_sites=s"   => \$flankingAttDir,
    "max_eval"      => \$maxEval,
    "force"         => \$force,
    "verbose"       => \$verbose,
    "help"          => \$help
    )
or die("Unknown option, try --help\n");

if ($help) {
    help();
    exit;
}

unless (-f $tablePath) {
    die "--dfam input is not a file: $!";
}

unless (-f $genomePath) {
    die "--genome input is not a file: $!";
}

unless (-d $tsvDir) {
    die "--tsv input is not a directory: $!";
}

unless (-d $chartDir) {
    die "--index_charts input is not a directory: $!";
}

unless (-d $flankingAttDir || not $flankingAttDir) {
    die "--att_sites input is not a directory: $!";
}

#extract bacteria name from table name
if ($tablePath =~ /([^\/]+)\.dfam/) {
    $dfamName = $1;
}
else {
    die "Can't extract .dfam name from table path $tablePath: $!";
}

print("dfam name: $dfamName\n");

#create directory that will contain index charts, unless it already exists
my $dir = "$chartDir/$dfamName";
print("new dir: $dir\n");
#if directory already exists, throw fatal error unless --force was specified
if (-e $dir && -d $dir) {
    unless ($force) {
        die "$dir already exists. Use --force to automatically overwrite it";
    }
    else {
        do_cmd("rm -r $dir");
        do_cmd("mkdir $dir");
    }
}
else {
    do_cmd("mkdir $dir");
}

$tsvPath = "$tsvDir/$dfamName.tsv";

parse_tables($tablePath, $genomePath, $tsvPath, $maxEval);

if ($verbose) {
    be_verbose();
}

sub parse_tables {
    #Uses a DFAM table corresponding to a genome to generate a .tsv file containing
    #entries for every prophage sequence detected, a file for any flanking att sites
    #found around those seqeunces, and nucleotide index charts that track which
    #nucleotides from the reference prophage strains were detected in the
    #bacterial genome.

    #Arguments:
        #0: path to DFAM table
        #1: Path to genome
        #2: path to output 'bed-like' file

    my $tablePath = $_[0];
    my $genomePath = $_[1];
    my $outputTSV = $_[2];
    my $maxEval = $_[3];
    my $count = 1;
    my %chartHash;
    my $chartPath;

    open(my $tableFile, "<", $tablePath) or die "Can't open $tablePath: $!";
    open(my $outputTSF, ">", $outputTSV) or die "Can't open $outputTSV: $!";

    while (my $line = <$tableFile>) {

        # Skip lines starting with #
        if ($line =~ /#/) {
        }

        #group 1: target name
        #group 2: match e-value
        #group 3: match start position on target (phage)
        #group 4: match end position on target (phage)
        #group 5: positive or negative strand
        #group 6: match start position on query (bacterial genome)
        #group 7: match end position on query (bacterial genome)
        #group 8: target model length (length of prophage)
        elsif ($line =~ /(.+?)\s+.+?\s+.+?\s+.+?\s+(.+?)\s+.+?\s+(.+?)\s+(.+?)\s+(.+?)\s+(.+?)\s+(.+?)\s+.+?\s+.+?\s+(\d+)\s/){
            my $name = $1;
            my $strand = $5;

            #Create a ViralSeq object. This will handle creating lines for the
            #output .tsv and detecting flanking att sites.
            my $seq = ViralSeq->new(
                name                => "$name",
                evalue              => $2,
                refSt               => $3,
                refEn               => $4,
                referenceSeqLength  => $8,
                strand              => $strand,
                gnSt                => $6,
                gnEn                => $7,
                flankingAttDir      => $flankingAttDir,
                genomePath          => $genomePath,
                verbose             => $verbose,
            );


            # if option was specified by user, run this
            if ($flankingAttDir) {
                $seq->findFlankingAtts();
            }

            if ($seq->evalue <= $maxEval) {
                if ($count == 1) {
                    print $outputTSF $seq->tableHeader() . "\n";
                }

                print $outputTSF $seq->tableLine() . "\n";

                $chartPath = "$chartDir/$dfamName/$name" . ".txt";

                #unless an entry for this prophage already exists in hash, generate
                #a new array for this sequence by putting each line into an array
                #index. Store this array in the hash
                unless (exists $chartHash{$name}) {
                    my @array = (0)x$seq->referenceSeqLength;
                    my $arrayRef = \@array;

                    #if file already exists, read in its contents and populate array
                    if (open(my $chartFile, "<", $chartPath)) {
                        # skip header
                        readline $chartFile;

                        my $index = 0;
                        while (my $line = <$chartFile>) {
                            chomp $line;
                            $arrayRef->[$index] = $line;
                            $index += 1;
                        }

                        close $chartFile;
                    }
                    $chartHash{$name} = $arrayRef;
                }

                #Each array index corresponds to a nucleotide in the reference strain.
                #Iterate over nucleotides present in genomic strain, incrementing
                #for each present in genome. To match array 0-indexing, nucleotide
                #indexes are offset by 1.
                #refSt = start location on reference seq, refEn = end location on reference seq
                for (my $i = $seq->refSt - 1; $i < $seq->refEn; $i++) {
                    my $size = $seq->referenceSeqLength;

                    if ($i < 0 || $i >= $size) {
                        die "Nucleotide index is outside of bounds 1-$size for $name";
                    }
                    my $arrayRef = $chartHash{$name};
                    $arrayRef->[$i] += 1;
                }
            #we only care about lines that matched the regex statement and had a good e-value, so count increments inside of if statements
            $count++;
            }
        }
    }

    close $tableFile;

    #since count is set to 1 initially for naming purposes, decrement by 1 for
    #reporting accuracy
    $count = $count - 1;

    my @hashKeys = keys %chartHash;
    print   "\nNumber of passing prophage matches in $tablePath: $count" . 
            "\n" . scalar @hashKeys . " strains had hits in $tablePath\n";

    #Print out index-based 'charts' where each index corresponds to a line
    foreach my $hashKey (@hashKeys) {
        my $chartDir = "$chartDir/$dfamName/$hashKey";
        my $chartPath = "$chartDir/$dfamName/$hashKey" . "Chart.txt";
        # create output dir
        mkdir $chartDir or die "Can't create output chart dir: $!\n";
        # use quotemeta to escape any non-word characters in path (pipe delimiters are fairly common, but don't play nice with unix)
        $chartPath = quotemeta($chartPath);
        open(my $chartOutput, ">", $chartPath) or die "Can't open $chartPath: $!";

        my @chartArray = @{$chartHash{$hashKey}};
        my $arraySize = @chartArray;
        my $chartLine = "$hashKey\n";

        foreach my $entry (@chartArray) {
            $chartLine .= "$entry\n";
        }

        print $chartOutput "$chartLine";
        close $chartOutput or die "Can't close $chartPath: $!";
    }
}

sub help {
    print "
#table_parser.pl: Extract information about viral sequences from DFAM tables

Basic options:
    --help: Displays this page
    --verbose: Print additional information about values held in
        variables and commands used by table_parser
    --force: Overwrite pre-existing directories in the index chart directory
        if they already exist. This will erase all files in these
        directories
    --max_eval: Maximum allowable match evalue for match to be used

Input options:
    --dfam: Path to input .dfam table
    --genome: Path to .fasta format bacterial genome corresponding to input
        .dfam file

Output options:
    --tsv: Path to .tsv directory. All .tsv file values are tab-delimited
    --index_charts: Path to nucleotide index chart directory. Each files in 
        this directory will be the length of its respective viral genome + 1
        (because the first line is a header line). Each line therefore
        corresponds to a position in the viral genome, with its int value
        being the number of times that position was found to have a match
        in the bacterial genome
    --attsites: Path to flanking att side directory. Att site detection is
        spotty at best in its current implementation- use of this feature
        is not currently recommended.

";
}

sub do_cmd {
    my $cmd = $_[0];

    if ($verbose) {
        print "Command: $cmd\n";
    }

    return `$cmd`;
}

sub be_verbose {
    my $path = `pwd`;

    print "\nTable: $tablePath\n";
    print "Genome: $genomePath\n";
    print "Bacteria Name: $dfamName\n";
    print "TSV Path: $tsvPath\n";
    print "Chart directory: $chartDir\n";
    print "Flanking att site directory: $flankingAttDir\n";
    print "Current directory: $path\n";
}
