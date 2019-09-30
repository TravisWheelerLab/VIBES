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
my $hmmPath = '';
my $tableDir = '';
my $scannedTableDir = '';
my $mvMatchless = '';
my $cpu = 0;
my $verbose = 0;
my $help = 0;

GetOptions (
    "hmm_db=s"          => \$hmmPath,
    "job_number=i"      => \$genomeNumber,
    "dfam_dir=s"        => \$tableDir,
    "scan_dir=s"        => \$scannedTableDir,
    "gen_dir=s"         => \$genomeDir,
    "mv_matchless=s"    => \$mvMatchless,
    "cpu=i"             => \$cpu,
    "verbose"           => \$verbose,
    "help"              => \$help
    )
or die("Unknown argument, try --help\n");

if($help) {
    help();
    exit
}

if ($cpu < 0) {
    die("--cpu must be at least 0\n")
}

#we assume every entry in the directory is a genome we're interested in
my @genomes = glob "$genomeDir/*";

#Extract file name from path, then run nhmmscan and dfamscan
my $genome = $genomes[$genomeNumber];
my @splitLine = split('/', $genome);
my $fileName;
my $tablePath;
my $scannedPath;

if ($splitLine[-1] =~ /(.+)\./) {
    $fileName = $1;
}
else {
    if ($verbose) {
        be_verbose();
    }
    die "Unable to extract file name from path: $genome\n";
}

$tablePath = "$tableDir/$fileName.dfam";
$scannedPath = "$scannedTableDir/$fileName.dfam";


my $stTime = [Time::HiRes::gettimeofday()];
do_cmd("nhmmscan --cpu $cpu --dfamtblout $tablePath $hmmPath $genome");
my $elapsed = Time::HiRes::tv_interval($stTime);

if (has_no_matches($tablePath) and $mvMatchless) {
    do_cmd("mv $tablePath $mvMatchless");
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
    print "Ref: $hmmPath\n";
    print "Table dir: $tableDir\n";
    print "Verbose: $verbose\n";
    print "Genome: $genome\n";
    print "Table Path: $tablePath\n";
    print "Scanned table path: $scannedPath\n";
    print "Time to complete nhmmscan: $elapsed seconds\n";
    print "Current directory: $path";
}

sub help {
    print "
#tableizer.pl: Use nhmmscan to search gen_dir genomes for hmm_db sequences or 
#genomes. Outputs results in DFAM tabular format and automatically scans with 
#dfamscan.pl, removing redundant hits. Intended to be run on a server cluster
#with a job manager assigning a number to each instance.
-------------------------------------------------------------------------------
--help: Prints this message

Required:
 --hmm_db <s>       Path to hmm db to search against genomes
 --job_number <i>   For n genomes, 0 <= i <= n-1
 --dfam_dir <s>     Path to directory of unscanned .dfam tables, which will
                    retain redundant hits
 --scan_dir <s>     Path to directory of scanned .dfam tables, which have 
                    redundant hits removed
 --gen_dir <s>      Path to directory of baterial genomes in .fasta format

Optional:
 --verbose          Prints debugging information
 --cpu <i>          Tells nhmmscan how many worker threads to use (min 0). Total
                    number of threads will be (1 + <i>)
 --mv_matchless <s> Moves any .dfam files without any entries to a folder
                    intended to hold matchless files.
";

}

sub do_cmd {
    my $cmd = $_[0];

    if ($verbose > 0) {
        print "$cmd\n";
    }
    my $res = `$cmd 2>&1`;

    if ($verbose > 0) {
        print "$res\n";
    }
    
    return ;
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