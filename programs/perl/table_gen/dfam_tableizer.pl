#!/usr/bin/env perl
#Takes in a path to a genome directory, a path to a file containing reference
#.hmms, paths to output directories, a suffix to append to the file name,
#and a boolean that sets the value of verbose.
#Returns a .dfam table for each genome in the directory.

#Thanks to Kaitlin Carey for help in understanding how to run things on the cluster!

use strict;
use warnings;
use Time::HiRes;
use Getopt::Long;

my $hmmPath = '';
my $genomePath = '';
my $tablePath = '';
my $scannedPath = '';
my $mvMatchless = '';
my $cpu = 0;
my $verbose = 0;
my $help = 0;

GetOptions (
    # Inputs
    "hmm_db=s"          => \$hmmPath,
    "genome=s"          => \$genomePath,

    # Outputs
    "dfam=s"            => \$tablePath,
    "scanned_dfam=s"    => \$scannedPath,

    # Options
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

unless(-f $hmmPath) {
    die "--hmm_db is not a file: $!";
}

if ($cpu < 0) {
    die("--cpu must be at least 0\n")
}

# make sure $genomePath exists before we go further
unless(-f $genomePath) {
    die "$genomePath is not a file: $!";
}

my $stTime = [Time::HiRes::gettimeofday()];
do_cmd("nhmmscan --cpu $cpu --dfamtblout $tablePath $hmmPath $genomePath");
my $elapsed = Time::HiRes::tv_interval($stTime);

if (has_no_matches($tablePath) and $mvMatchless) {
    do_cmd("mv $tablePath $mvMatchless");
}
else {
    use File::Basename;
    my $dirname = dirname(__FILE__);
    do_cmd("perl $dirname/dfamscan.pl --dfam_infile  $tablePath --dfam_outfile $scannedPath");
}

if ($verbose) {
    be_verbose();
}

sub be_verbose {
    my $path = `pwd`;

    print "Ref: $hmmPath\n";
    print "Genome: $genomePath\n";
    print "Dfam: $tablePath\n";
    print "Scanned Dfam: $scannedPath\n";
    print "Verbose: $verbose\n";
    print "Time to complete nhmmscan: $elapsed seconds\n";
    print "Current directory: $path";
}

sub help {
    print "
# $0:
#Use nhmmscan to search input genome for hmm_db sequences/genomes. 
#Outputs results in DFAM tabular format and automatically scans with 
#dfamscan.pl, removing redundant hits. Intended to be run on a server cluster
#with a job manager assigning a number to each instance.
-------------------------------------------------------------------------------
--help: Prints this message

Required:
 --hmm_db <s>       Path to hmm db to search against genomes
 --dfam_dir <s>     Path to directory of unscanned .dfam tables, which will
                    retain redundant hits
 --scan_dir <s>     Path to directory of scanned .dfam tables, which have 
                    redundant hits removed
 --genome   <s>      Path to bacterial genome to search against

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

    # catch error thrown by backticks command
    if ($? != 0) {
        die "command `$cmd` returned error: $res";
    }

    if ($verbose > 0) {
        print "$res\n";
    }
    
    return $res;
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
