#!/usr/bin/env perl

use strict;
use warnings;
use Time::HiRes;
use Getopt::Long;

#command line argument variables that can be provided by user
my $genomeNumber;
my $genomeDir = '';
my $referenceHMMPath = '';
my $outputDir  = '';
my $verbose = 0;
my $help = 0;
my $force = 0;

#other variables
my $outputPath = '';
my $genome = '';
my $elapsed = '';

GetOptions (
    "hmm=s"         => \$referenceHMMPath,
    "outputdir=s"   => \$outputDir,
    "jobnumber=i"   => \$genomeNumber,
    "genomedir=s"   => \$genomeDir,
    "verbose"       => \$verbose,
    "help"          => \$help,
    "force"         => \$force
    )
or die("Unknown argument, try --help\n");

if ($help) {
    help();
}
else {
    run_phmmer();

    if($verbose) {
        be_verbose();
    }
}

sub run_phmmer {
    #put all files in directory into an array
    my @genomes = glob "$genomeDir/*";
    $genome = $genomes[$genomeNumber];

    #Extract file name from genome path
    my @splitLine = split('/', $genome);
    my $fileName;

    if ($splitLine[-1] =~ /(.+?)\./) {
        $fileName = $1;
    }
    else {
        be_verbose();
        die "Unable to extract file name from path: $genome\n";
    }

    $outputPath = "$outputDir/$fileName.tbl";

    #Ensure output file doesn't already exist, and that user wants to overwrite it
    #if it does. Time how long it takes PHMMERT to run, which can be reported with
    # --verbose
    checkFile($outputPath);
    my $stTime = [Time::HiRes::gettimeofday()];
    do_cmd("phmmert --tblout $outputPath $referenceHMMPath $genome");
    $elapsed = Time::HiRes::tv_interval($stTime);
}

# execute command at command line, printing command if --verbose
sub do_cmd {
    my $cmd = $_[0];

    if ($verbose > 0) {
        print "$cmd\n";
    }

    return `$cmd`;
}

sub help { 
    print "
# cluster_phmmert.pl: runs phmmert with --tblout on genomes

Usage: perl cluster_phmmert [options] --hmm [path] --genomedir [path]
    --jobnumber [int] --outputdir [path]

Options:
    --help: Displays this help page
    --verbose: Prints any commands run and other info. Intended for debugging
    --force: With this selected, any files with the same name as an output file
        in the output directory will be automatically overwritten.

Input:
    --genomedir: Path to directory populated with genomes to run phmmert against
    --jobnumber: Int assigned by server cluster, which specifies a genome to use
    --hmm: Path to reference HMM to use as phmmert query

Output:
    --outputdir: Path to directory to be filled with output phmmert tables"
}

# Check if a file already exists. If yes, and --force hasn't been specified by the user,
# print an error to a log file and STDERR, then stop the program. If --force has
# been specified, then overwrite the file
sub checkFile {
   my $file = $_[0];

   if (-f $file && !$force) {
        open(my $errorlog, '>>', "cluster_phmmert_errors.txt") or die "Could not open file '$file' $!";
        print $errorlog "$file already exists! To overwrite any files that already exist, rerun with --force option.\n\n";
        die "$file already exists! To overwrite any files that already exist, rerun with --force option.\n";
  }
}

sub be_verbose {
    my $path = `pwd`;

    print "Number: $genomeNumber\n";
    print "GenomeDir: $genomeDir\n";
    print "RefHMM: $referenceHMMPath\n";
    print "Output dir: $outputDir\n";
    print "Genome: $genome\n";
    print "Time to complete phmmert: $elapsed seconds\n";
    print "Current directory: $path";
    print "Verbose: $verbose\n";
    print "Force: $force\n";
}