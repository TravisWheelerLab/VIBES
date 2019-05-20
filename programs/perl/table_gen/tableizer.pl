#!/usr/bin/env perl
#Takes in a path to a genome directory, a path to a file containing reference
#.hmms, paths to output directories, a suffix to append to the file name,
#and a boolean that sets the value of verbose.
#Returns a .dfam table for each genome in the directory.

#Credit to Kaitlin Carey for help in understanding how to run things on the cluster!

use strict;
use warnings;
use Time::HiRes;
use Getopt::Long;

my $genomeNumber;
my $genomeDir = '';
my $referencePath = '';
my $tableDir = '';
my $scannedTableDir = '';
my $suffix = '';
my $rmMatchless = '';
my $verbose = 0;
my $help = 0;

GetOptions (
    "hmm_db=s"          => \$referencePath,
    "job_number=i"      => \$genomeNumber,
    "dfam_dir=s"        => \$tableDir,
    "output_dir=s"      => \$scannedTableDir,
    "bac_dir=s"         => \$genomeDir,
    "suffix=s"          => \$suffix,
    "rm_matchless=s"      => \$rmMatchless,
    "verbose"           => \$verbose,
    "help"              => \$help
    )
or die("Unknown argument, try --help\n");

if($help) {
    help();
    exit
}

#we assume every entry in the directory is a genome we're interested in
my @genomes = glob "$genomeDir/*";

#Extract file name from path, then run nhmmscan and dfamscan
my $genome = $genomes[$genomeNumber];
my @splitLine = split('/', $genome);
my $fileName;

if ($splitLine[-1] =~ /(.+)\./) {
    $fileName = $1;
}
else {
    if ($verbose) {
        be_verbose();
    }
    die "Unable to extract file name from path: $genome\n";
}

my $tablePath = "$tableDir/$fileName" . "_$suffix.dfam";
my $scannedPath = "$scannedTableDir/$fileName" . "_$suffix" . "_scanned.dfam";

my $stTime = [Time::HiRes::gettimeofday()];
do_cmd("nhmmscan --cpu 1 --dfamtblout $tablePath $referencePath $genome");
my $elapsed = Time::HiRes::tv_interval($stTime);

if (has_no_matches($tablePath) and $rmMatchless) {
    do_cmd("mv $tablePath $rmMatchless");
}
else {
    do_cmd("perl dfamscan.pl --dfam_infile  $tablePath --dfam_outfile $scannedPath");
}

if ($verbose) {
    be_verbose();
}

sub be_verbose {
    my $path = `pwd`;

    print "Number: $genomeNumber\n";
    print "GenDir: $genomeDir\n";
    print "Ref: $referencePath\n";
    print "Table dir: $tableDir\n";
    print "Suffix: $suffix\n";
    print "Verbose: $verbose\n";
    print "Genome: $genome\n";
    print "Table Path: $tablePath\n";
    print "Scanned table path: $scannedPath\n";
    print "Time to complete nhmmscan: $elapsed seconds\n";
    print "Current directory: $path";
}

sub help {
    print "
#tableizer.pl: Convert genomes into DFAM format tabular output. Intended to be
run on a server cluster with a job manager assigning a number to each instance.
-------------------------------------------------------------------------------
--help: Prints this message

Required:
 --hmm_db <s>       Path to hmm db to search against genomes
 --job_number <i>   For n genomes, 0 <= i <= n-1
 --dfam_dir <s>     Path to directory of unscanned .dfam tables, which will
                    retain redundant hits
 --output_dir <s>   Path to directory of scanned .dfam tables, which have 
                    redundant hits removed
 --bac_dir <s>      Path to directory of baterial genomes in .fasta format
 --suffix <s>       Suffix to append to file names
                    (i.e. strain_name_suffix.dfam)

Optional:
 --verbose          Prints debugging information
";

}

sub do_cmd {
    my $cmd = $_[0];

    if ($verbose > 0) {
        print "$cmd\n";
    }

    return `$cmd`;
}

# checks if a .dfam file has no matches. Returns 1 if no matches are found, 0 if matches are present
sub has_no_matches {
    my $file = $_[0];
    my $isMatchless = 1;

    # iterate through lines in file, looking for any lines that don't start with #.
    # # is used to signify lines that don't contain data
    open (my $fileHandle, "<", $file);

    while (my $line = <$fileHandle>) {
        # if the line doesn't start with #, the file contains at least 1 match
        if (not $line =~ /^#/) {
            $isMatchless = 0;
        }
    }

    return $isMatchless;
}
