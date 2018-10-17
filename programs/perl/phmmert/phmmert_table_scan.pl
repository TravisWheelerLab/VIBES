#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use List::Util qw( min max );
#credit to user ikegami at https://stackoverflow.com/questions/10701210/how-to-find-maximum-and-minimum-value-in-an-array-of-integers-in-perl
#for suggesting the use of List::Util

#command line argument variables that can be set by user
my $tableDir = '';
my $output = '';
my $force = 0;
my $rmEmpty = 0; #if a table file is found to be empty, this will signal whether or not to delete it
my $help = 0;
my $verbose = 0;

#other variables
my %bestHitHash = ();
my $headerLine = "genome,best match e-value,best match length";

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
            my $minVal = "inf"; #since we want the minimum e-value, we set the starting
            #value to be the largest possible value
            my $bestLength = 0;

            #go through lines in the table, keeping the length and e-value of the best match
            while (my $line = <$tableHandle>) {
                #if the line matches a "hit" line, grab the e-value
                if ($line =~ /.+?(\d+)\s+(\d+)\s+\d+\s+\d+\s+.+?\s+(.+?)\s+.+?\s+.+?$/) {
                    my $alignStart = $1;
                    my $alignEnd = $2;
                    my $eValue = $3;

                    my $length = abs($alignStart - $alignEnd); #length could be negative due to
                    #negative strand

                    if ($eValue < $minVal) {
                        $minVal = $eValue;
                        $bestLength = $length;
                    }
                }

            }

            #since no e-value should be infinity, we know a line was successfully read in,
            #and that the file was therefore not empty
            if ($minVal != "inf") {
                my @values = [];
                my $valuesRef = \@values;

                $values[0] = $minVal;
                $values[1] = $bestLength;

                $bestHitHash{$genome} = $valuesRef;
            }
            elsif ($rmEmpty) {
                do_cmd("rm $table");
            }
        }
    }
}

#Prints best contents of non-empty tables to user-specified output file
sub print_results {
    open(my $outputHandle, ">", $output) or die "Can't open output file $output: $!";
    print $outputHandle "$headerLine\n";

#hash keys are genome names, so sort and iterate through them to print to output
    foreach my $genome (sort keys %bestHitHash) {
        my $valuesRef = $bestHitHash{$genome};


        print $outputHandle "$genome," . $valuesRef->[0] . "," . $valuesRef->[1] . "\n";
    }
}