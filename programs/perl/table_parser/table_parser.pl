#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib $FindBin::Bin;
# import classes
use ViralSeq;
use Getopt::Long;

#default values reflect file structure I set up on the cluster
my $refProphageDir = "../prophage_seq";
my $tableDir = "../tables/scanned";
my $genomeDir = "../ps_genomes";
my $tsvDir = "../tsv";
my $chartDir = "../nuc_chart";
my $flankingAttDir = "../flanking_att_sites";
my $jobNumber = 0;
my $verbose = ''; #default false value
my $help = ''; #^^
my $force = 0; #^^
my $genomeName;
my $genomePath;
my $tsvPath;

GetOptions (
    "prophage=s"    => \$refProphageDir,
    "tables=s"      => \$tableDir,
    "genomes=s"     => \$genomeDir,
    "tsv=s"         => \$tsvDir,
    "indexcharts=s" => \$chartDir,
    "attsites=s"    => \$flankingAttDir,
    "jobnumber=i"   => \$jobNumber,
    "force"         => \$force,
    "verbose"       => \$verbose,
    "help"          => \$help
    )
or die("Unknown option, try --help\n");

if ($help) {
    help();
    exit;
}

my @tables = glob "$tableDir/*" or die "Can't find $tableDir: $!";
my $table = $tables[$jobNumber];

#extract genome name from table name
if ($table =~ /([^\/]+)_viral_scanned\.dfam/) {
    $genomeName = $1;
}
else {
    die "Can't extract genome name from table path $table: $!";
}

#create directory that will contain index charts, unless it already exists
my $dir = "$chartDir/$genomeName";
#if directory already exists, throw fatal error unless --force was used
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

$tsvPath = "$tsvDir/$genomeName.tsv";
$genomePath = "$genomeDir/$genomeName.fna";

parse_tables($table, $genomePath, $tsvPath);

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
    my $count = 1;
    my %chartHash;
    my $chartPath;

    open(my $tableFile, "<", $tablePath) or die "Can't open $tablePath: $!";
    open(my $outputTSF, ">", $outputTSV) or die "Can't open $outputTSV: $!";

    while (my $line = <$tableFile>) {

        # Skip lines starting with #
        if ($line =~ /#/) {
        }

        #Well folks, it's already time for a monsterous, unreadable Regex line.
        #To summarize its twisted inner workings: Group 1 captures the name of the
        #strain, group 2 nabs the start index of the reference strain hit,
        #group 3 grabs the end index of the reference strain hit, group 4 obtains
        #the + or - sign indicating which strand the sequence is one, group 5 nicks
        #the start index in the bacterial genome, and group 6 contains the end index
        #in the bacterial genome.
        elsif ($line =~ /(.+?)\s+.+?\s+.+?\s+.+?\s+.+?\s+.+?\s+(.+?)\s+(.+?)\s+(.+?)\s+(.+?)\s+(.+?)\s/){
            my $name = $1;
            my $strand = $4;
            my $isPos = 0;

            if ($strand eq "+"){
                $isPos = 1;
            }
            #Create a ViralSeq object. This will handle creating lines for the
            #output .tsv and detecting flanking att sites.
            my $seq = ViralSeq->new(
                name                => "$name.$5",
                refSt               => $2,
                refEn               => $3,
                referenceSeqPath    => "$refProphageDir/$name.fasta",
                isPos               => $isPos,
                gnSt                => $5,
                gnEn                => $6,
                flankingAttDir      => $flankingAttDir,
                genomePath          => $genomePath,
                verbose             => $verbose,
            );

            $seq->findFlankingAtts();

            if ($count == 1) {
                print $outputTSF $seq->tableHeader() . "\n";
            }

            print $outputTSF $seq->tableLine() . "\n";

            $chartPath = "$chartDir/$genomeName/$name" . "Chart.txt";

            #unless an entry for this prophage aleady exists in hash, generate
            #a new array for this sequence by putting each line into an array
            #index. Store this array in the hash
            die "Change how new arrays are generated- since we're creating a new
            dir for each genome, if an entry doesn't exist, then we ";
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
            #we only care about lines that matched the regex statement, so count
            #increments inside of if statement
            $count++;
        }
    }

    close $tableFile;

    #since count is set to 1 initially for naming purposes, decrement by 1 for
    #reporting accuracy
    $count = $count - 1;

    my @hashKeys = keys %chartHash;
    print   "\nNumber of prophage sequences in $tablePath: $count" . 
            "\nOf 50 reference sequences, " . scalar @hashKeys . " had hits in $tablePath\n";

    #Print out index-based 'charts' where each index corresponds to a line
    foreach my $hashKey (@hashKeys) {
        $chartPath = "$chartDir/$genomeName/$hashKey" . "Chart.txt";
        open(my $chartOutput, ">", $chartPath) or die "Can't open $chartPath: $!";

        my @chartArray = @{$chartHash{$hashKey}};
        my $chartLine = "$hashKey\n";

        foreach my $entry (@chartArray) {
            $chartLine .= "$entry\n";
        }

        print $chartOutput "$chartLine";
        close $chartOutput;
    }
}

sub help {
    print "
    #table_parser.pl: Extract information about viral sequences from DFAM tables
    #Input and output options have default values, but can be specified by user

    Basic options:
        --help: Displays this page
        --verbose: Print additional information about values held in
            variables and commands used by table_parser
        --force: Overwrite pre-existing directories in the index chart directory
            if necessary

    Input options:
        --prophage: Path to directory containing reference prophage
            strains
        --tables: Path to DFAM table directory
        --genomes: Path to directory with bacterial genomes
        --jobnumber: Integer provided by the cluster job manager that tells
            table_parser which table it should run on

    Output options:
        --tsv: Path to .tsv directory
        --indexcharts: Path to nucleotide index chart directory
        --attsites: Path to flanking att side directory

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

    print "\nTable directoy: $tableDir\n";
    print "Reference prophage sequence directory: $refProphageDir\n";
    print "Genome directory: $genomeDir\n";
    print "Genome path: $genomePath\n";
    print "Genome Name: $genomeName\n";
    print "TSV Path: $tsvPath\n";
    print "Chart directory: $chartDir\n";
    print "Flanking att site directory: $flankingAttDir\n";
    print "Current directory: $path\n";
}
