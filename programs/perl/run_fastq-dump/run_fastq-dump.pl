#!/usr/bin/env perl

# takes in a file of Sequence Read Archive (SRA) accession numbers and an output folder. Runs fastq-dump on each acc number, placing resultant
# .fastq file in output folder
use strict;
use warnings;

use Getopt::Long;

# Command-line options
my $inputPath = '';
my $outputDir = '';
my $logPath = '';
my $help = 0;
my $verbose = 0;
my $force = 0;

GetOptions (
    "input=s"       => \$inputPath,
    "output=s"      => \$outputDir,
    "log=s"         => \$logPath,
    "help"          => \$help,
    "force"         => \$force,
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
        print $log "=====\n";
        my $date = do_cmd("date");
        print $log "$date\n";
        print $log "
Input path: $inputPath
Output dir: $outputDir
Help: $help
Verbose: $verbose
Force: $force
";
    }

    while (my $accession = <$accessionFile>) {
        # strip whitespace characters such as \n from accession number. 'g' flag instructs regex to keep going after first match
        $accession =~ s/\s//g;

        # unless an entry already exists for accession and --force wasn't enabled
        unless(-e "$outputDir/$accession.fastq" && !$force) {
            do_cmd("fastq-dump -O $outputDir $accession", $log);
        }
        else {
            my $string = "File already found for $accession\n";
            print($string);
            if ($logPath) {
              print $log $string;
            }
        }
    }
    close $accessionFile;
    close $log;
}

sub do_cmd {
    my $cmd = $_[0];
    my $log= $_[1];

    if ($verbose > 0) {
        print "$cmd\n";
    }
    my $res = `$cmd 2>&1`;

    if ($verbose > 0) {
        print "$res\n";
    }

    # If logPath set by user
    if ($log) {
        print $log "$cmd\n";
        print $log "$res\n\n";
    }

    
    return $res;
}

sub help {
    print "#$0
#Reads in SRA accession IDs from file and downloads associated .fastq files with
#fastq-dump
--------------------------------------------------------------------------------
Input:
    --input <s>     Path to input .txt files of SRA accession IDs, separated by 
                    newline characters
    
Output:
    --output <s>    Path to output directory where .fastq files will be stored.
                    Each output .fastq file will be named after its accession #

    --log <s>       Path to log file, where commands run by program and
                    fastq-dump commandline output are stored

Misc:
    --help          Display this help pane
    --force         Automatically overwrite output .fastq files, if they already
                    exist
    --verbose       Print all commands run to commandline output
"
}