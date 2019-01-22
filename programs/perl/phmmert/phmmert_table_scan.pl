#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

#command line argument variables that can be set by user
my $tableDir = '';
my $output = '';
my $force = 0;
my $rmEmpty = 0; #if a table file is found to be empty, this will signal whether or not to delete it
my $help = 0;
my $verbose = 0;

#other variables
my %bestHitHash = ();
my $headerLine = "Genome,Best match e-value,Best match length";
my $minVal = "inf"; #this will act as the minimum acceptable e-value

GetOptions (
    "tabledir=s"    => \$tableDir,
    "output=s"      => \$output,
    "help"          => \$help,
    "rmempty"       => \$rmEmpty,
    "force"         => \$force,
    "verbose"       => \$verbose
    )
or die("Unknown argument, try --help\n");

if ($help) {
    help();
}
else {
    scan_files();
    print_results();
}

#Scans through all files in the table directory, which are expected to be
#phmmert tables. After scanning through a file, we determine whether or not
#it's empty and delete it if the user specifies to do so. Non-empty files 
sub scan_files {
    my @tables = glob "$tableDir/*";
    my $genome = '';

    #Grab genome name from file path
    foreach my $table (@tables) {
        if ($table =~ /([^\/]+)\.tbl/) {
            $genome = $1;
        }
        else {
            die "Unable to parse genome name from file $table\n";
        }

        #double-check that "table" is indeed a file, then read in contents
        if (-f $table) {
            open(my $tableHandle, "<", $table) or die "Can't open $table: $!";

            my $numberHits = 0; #we track the number of hits detected, to determine
            #whether a table is empty
            my $bestMatch = $minVal; #set best match to minimum allowable e-value
            my $isEmpty = 1; #used to track whether a file contains any entries at all
            my $hasAcceptableEntry = 0; #boolean to track whether or not we find an entry with a good enough e-value
            my $bestLength = 0;

            #go through lines in the table, keeping the length and e-value of the best match
            while (my $line = <$tableHandle>) {
                #if the line matches a "hit" line, grab the e-value
                if ($line =~ /.+?(\d+)\s+(\d+)\s+\d+\s+\d+\s+.+?\s+(.+?)\s+.+?\s+.+?$/) {
                    my $alignStart = $1;
                    my $alignEnd = $2;
                    my $eValue = $3;

                    $isEmpty = 0;
                    my $length = abs($alignStart - $alignEnd) + 1; #length could be negative due to
                    #negative strand

                    if ($eValue <= $bestMatch) {
                        $bestMatch = $eValue;
                        $bestLength = $length;
                        $hasAcceptableEntry = 1;
                    }
                }

            }

            #save best e-value (if any) in hash
            if ($hasAcceptableEntry) {
                my @values = [];
                my $valuesRef = \@values;

                $values[0] = $bestMatch;
                $values[1] = $bestLength;

                $bestHitHash{$genome} = $valuesRef;
            }
            elsif ($rmEmpty && $isEmpty) {
                do_cmd("rm $table");
            }
        }
    }
}

sub help {
    print "\n
#phmmert_table_scan.pl: Scans through a directory of phmmert tables, then prints
#the e-value and length of the best match of each table to a file in .csv format
#Expects all files in the input directory to be phmmert table files

Usage: perl phmmert_table_scan.pl [options] --tabledir [path] --output [path]

Options:
    --help: Displays this help page
    --verbose: Prints out console commands used by program
    --rmempty: Automatically delete tables that don't have any entries
    --force: Overwrite pre-existing files instead of stopping the program

Input:
    --tabledir: Path to a directory containing phmmert table files. All files in
        this directory are expected to be phmmert tables.

Output:
    --output: Path to desired location of output .csv file
    \n";
}

# Check if a file already exists. If yes, and --force hasn't been specified by the user,
# print an error to a log file and STDERR, then stop the program. If --force has
# been specified, then overwrite the file
sub checkFile {
   my $file = $_[0];

   if (-f $file && !$force) {
        open(my $errorlog, '>>', "phmmert_table_scan_errors.txt") or die "Could not open file '$file' $!";
        print $errorlog "$file already exists! To overwrite any files that already exist, rerun with --force option.\n\n";
        die "$file already exists! To overwrite any files that already exist, rerun with --force option.\n";
  }
}

#Prints best contents of non-empty tables to user-specified output file
sub print_results {
    checkFile($output);
    open(my $outputHandle, ">", $output) or die "Can't open output file $output: $!";
    print $outputHandle "$headerLine\n";

#hash keys are genome names, so sort and iterate through them to print to output
    foreach my $genome (sort keys %bestHitHash) {
        my $valuesRef = $bestHitHash{$genome};


        print $outputHandle "$genome," . $valuesRef->[0] . "," . $valuesRef->[1] . "\n";
    }
}

# execute command at command line, printing command if --verbose
sub do_cmd {
    my $cmd = $_[0];

    if ($verbose > 0) {
        print "$cmd\n";
    }

    return `$cmd`;
}