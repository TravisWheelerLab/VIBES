#!/usr/bin/env perl

# takes in a file of Sequence Read Archive (SRA) accession numbers and an output folder. Runs fastq-dump on each acc number, placing resultant
# .fastq file in output folder
use strict;
use warnings;

use Getopt::Long;

# Command-line options
my $inputPath = '';
my $outputDir = '';
my $jobs = 0;
my $logPath = '';
my $help = 0;
my $verbose = 0;

GetOptions (
    "input=s"       => \$inputPath,
    "output=s"      => \$outputDir,
    "jobs=i"        => \$jobs,
    "log=s"         => \$logPath,
    "help"          => \$help,
    "verbose"       => \$verbose
    )
or die("Unknown option, try --help\n");

if ($help) {
    help();
    exit;
}
else {
    open my $accessionFile, "<$inputPath";
    my $log;
    if ($logPath) {
        open $log, ">>$logPath";
        my $date = do_cmd("date");
        logPrint("=====
$date
Input path: $inputPath
Output dir: $outputDir
Help: $help
Verbose: $verbose

", $log);
    }

    # use GNU parallel to run fastq-dump on every accID in the input .txt file,
    # saving resultant .fastqs to $outputDir
    do_cmd("cat $inputPath | parallel -j $jobs fastq-dump  --gzip -O $outputDir {.}", $log);

    if ($logPath) {
       close $log;
    }
}

sub do_cmd {
    my $cmd = $_[0];
    my $log= $_[1];

    # print command to log before execution, ensuring it's saved if program is
    # interrupted. This is particularly important because we want it to be
    # apparent if a .fastq file didn't finish downloading, in which case its
    # command won't be followed by a report of reads downloaded in the log
    if ($verbose > 0) {
        print "$cmd\n";
    }
    logPrint("$cmd\n", $log);

    my $res = `$cmd`;

    if ($verbose > 0) {
        print "$res\n";
    }

    logPrint("$res\n", $log);

    
    return $res;
}

sub help {
    print "#$0
#Uses GNU parallel to read accession IDs from input file, run multiple instances
#of fastq-dump, and save resultant .fastq files in output dir. Can automatically
#compress output dir with tar once all fastq-dump runs complete
--------------------------------------------------------------------------------
Input:
    --input <s>     Path to input .txt files of SRA accession IDs, separated by 
                    newline characters
    --jobs <i>      Number of instances of fastq-dump to be run simultaneously.
    
Output:
    --output <s>    Path to output directory where .fastq files will be stored.
                    Each output .fastq file will be named after its accession #

    --log <s>       Path to log file, where commands run by program and
                    fastq-dump commandline output are store

Misc:
    --help          Display this help pane
    --verbose       Print all commands run to commandline output
";
}

# Print string to log
sub logPrint {
    my $string = $_[0];
    my $log = $_[1];

    if ($log) {
        print $log "$string";
    }
}