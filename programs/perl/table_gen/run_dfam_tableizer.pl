#!/usr/bin/env perl
#Runs dfam_tableizer.pl using GNU Parallel instead of a job manager
use strict;
use warnings;
use Time::HiRes;
use Getopt::Long;

my $jobs;
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
    "jobs=i"            => \$jobs,
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
    die("--cpu must be at least 0\n");
}

#we assume every entry in the directory is a genome we're interested in
my @genomes = glob "$genomeDir/*";
my $seqLast = @genomes - 1;

my $stTime = [Time::HiRes::gettimeofday()];
do_cmd("seq 0 $seqLast | parallel -j $jobs perl dfam_tableizer.pl --hmm_db $hmmPath --job_number {.} --dfam_dir $tableDir --scan_dir $scannedTableDir --gen_dir $genomeDir --cpu $cpu --verbose $verbose");
my $elapsed = Time::HiRes::tv_interval($stTime);

if ($verbose) {
    be_verbose();
}

sub be_verbose {
    my $path = `pwd`;

    print "GenDir: $genomeDir\n";
    print "Ref: $hmmPath\n";
    print "Table dir: $tableDir\n";
    print "Verbose: $verbose\n";
    print "Table dir: $tableDir\n";
    print "Scanned table dir: $scannedTableDir\n";
    print "Time to complete nhmmscan: $elapsed seconds\n";
    print "Current directory: $path";
}

sub help {
    print "
# $0:
# In cases where there's not a job manager to handle running multiple instances
# of dfam_tableizer simultanrously, this program will use GNU Parallel to do so
--------------------------------------------------------------------------------
--help: Prints this message

Required:
 --hmm_db <s>       Path to hmm db to search against genomes
 --jobs <i>         Maximum instances of dfam_tableizer we want to run
                    simultaneously
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
    
    return $res;
}
