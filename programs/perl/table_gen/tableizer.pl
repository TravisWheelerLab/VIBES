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

#Genome number from server starts at one, so decrement to match array indexing
my $genomeNumber;
my $genomeDir = '';
my $referencePath = '';
my $tableDir = '';
my $scannedTableDir = '';
my $suffix = '';
my $verbose = '';
my $help = 0;

GetOptions (
    "prophageinput=s"   => \$referencePath,
    "jobnumber=i"         => \$$genomeNumber,
    "tableinput=s"      => \$tableDir,
    "outputdirectory=s" => \$scannedTableDir,
    "genomes=s"         => \$genomeDir,
    "suffix=s"          => \$suffix,
    "verbose"           => \$verbose,
    "help"              => \$help
    )
or die("Unknown argument, try --help\n");

die "Implement help()!\n";

#we assume every entry in the directory is a genome we're interested in
my @genomes = glob "$genomeDir/*";

#Extract file name from path, then run nhmmscan and dfamscan
my $genome = $genomes[$genomeNumber];
my @splitLine = split('/', $genome);
my $fileName;

if ($splitLine[-1] =~ /(.+?)\./) {
    $fileName = $1;
}
else {
    be_verbose();
    die "Unable to extract file name from path: $genome\n";
}

my $tablePath = "$tableDir/$fileName" . "_$suffix.dfam";
my $scannedPath = "$scannedTableDir/$fileName" . "_$suffix" . "_scanned.dfam";

my $stTime = [Time::HiRes::gettimeofday()];
do_cmd("nhmmscan --cpu 1 --dfamtblout $tablePath $referencePath $genome");
my $elapsed = Time::HiRes::tv_interval($stTime);
do_cmd("perl dfamscan.pl --dfam_infile  $tablePath --dfam_outfile $scannedPath");

if ($verbose > 0) {
    be_verbose();
}

sub be_verbose {
    my $path = `pwd`;

    print "Number: $genomeNumber\n";
    print STDERR "Number: $genomeNumber\n";
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
    #tableizer.pl: Convert genomes into DFAM format tabular output";

}

sub do_cmd {
    my $cmd = $_[0];

    if ($verbose > 0) {
        print "$cmd\n";
    }

    return `$cmd`;
}
