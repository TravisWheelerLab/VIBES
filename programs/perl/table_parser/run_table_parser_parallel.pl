#!/usr/bin/env perl
#Runs table_parser.pl using GNU Parallel instead of a job manager
use strict;
use warnings;
use Time::HiRes;
use Getopt::Long;

my $refProphageDir = "";
my $tableDir = "";
my $genomeDir = "";
my $tsvDir = "";
my $chartDir = "";
my $flankingAttDir = "";
my $jobNumber = 0;
my $suffix = '';
my $verbose = ''; #default false value
my $help = ''; #^^
my $force = 0; #^^
my $genomeName;
my $genomePath;
my $tsvPath;
my $maxEval = 1e-5;
my $jobs = 2;

GetOptions (
    "prophage=s"    => \$refProphageDir,
    "dfam=s"        => \$tableDir,
    "bac_genomes=s" => \$genomeDir,
    "tsv=s"         => \$tsvDir,
    "index_charts=s"=> \$chartDir,
    "att_sites=s"   => \$flankingAttDir,
    "jobs=i"        => \$jobs,
    "suffix=s"      => \$suffix,
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

#we assume every entry in the directory is a genome we're interested in
my @tables = glob "tableDir/*";
# should give us length of @tables minus 1 (zero-indexed last index)
my $seqLast = @tables - 1;

# set first part of command string, which should always be the same: we use seq to generate all numbers from 0 to len(@tables) - 1, 
my $cmdString = "seq 0 $seqLast | parallel --dryrun -j $jobs perl table_parser.pl";

# depending on set flags, concatenate same flags to cmdString
if (defined $verbose) {
    $cmdString .= " --verbose";
}

if (defined $force) {
    $cmdString .= " --force";
}

if (defined $suffix) {
    $cmdString .= "--suffix $suffix";
}

#perl table_parser.pl --verbose --force --jobNumber $MOAB_JOBARRAYINDEX --prophage ../seq/prophage_seq/ --tables ../tables/dfam/2019-06-14_complete_scanned/ 
#--genomes ../seq/pseudomonas/ --tsv ../tables/tsv/2019-06-19_complete --indexcharts ../nuc_chart/2019-06-19_complete --suffix _scanned
$cmdString .= "--prophage $refProphageDir --tables $tableDir --bac_genomes $genomeDir --tsv $tsvDir --index_charts";

do_cmd($cmdString);

my $stTime = [Time::HiRes::gettimeofday()];
do_cmd(" -");
my $elapsed = Time::HiRes::tv_interval($stTime);

if ($verbose) {
    be_verbose();
}

sub be_verbose {
    my $path = `pwd`;

    print "GenDir: $genomeDir\n";
    print "Table dir: $tableDir\n";
    print "Verbose: $verbose\n";
    print "Table dir: $tableDir\n";
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

    if ($verbose) {
        print "$cmd\n";
    }
    my $res = `$cmd 2>&1`;

    if ($verbose) {
        print "$res\n";
    }
    
    return $res;
}
