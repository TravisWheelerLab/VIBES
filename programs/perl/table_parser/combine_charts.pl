#!/usr/bin/env perl

use strict;
use warnings;
use File::Find;
use Getopt::Long;

#Command-line options
my $inputDir = '';
my $outputDir = '';
my $help = '';
my $verbose = '';
my $force = '';

#internal variables
my %phageHash;

GetOptions (
    "inputdir=s"    => \$inputDir,
    "outputdir=s"   => \$outputDir,
    "help"          => \$help,
    "verbose"       => \$verbose,
    "force"         => \$force
    )
or die("Unknown option, try --help");

if ($help) {
    help();
}
else {
    #iterate through all files and directories in dir, giving them as input to scan_index_file()
    find(\&scan_index_file, $inputDir);

    #test printing statement
    # for my $key (keys %phageHash) {
    #     print "Key: $key\n";
    #     print "@{$phageHash{$key}}\n";
    # }

    print_combined_charts();
}

#use File::Find to iterate through the directory of genome directories, adding up the results for all of the phage index counts
sub scan_index_file {
    my $prophage = '';
    my $file = $_;

    #unless the 'file' we're examining is a direc1tory, 
    unless (-d) {
        #grab the prophage name from file name
        if ($file =~ /(.+?)Chart.txt/) {
            $prophage = $1;
        }
        else {
            die "Unable to parse file name: $file\n";
        }

        my $arrayRef = undef;

        #if the hash doesn't have an entry for this prophage, create one
        unless (exists $phageHash{$prophage}) {
            my @array = ();
            $arrayRef = \@array;
            $phageHash{$prophage} = $arrayRef;
        }
        #else use pre-existing hash entry
        else {
            $arrayRef = $phageHash{$prophage};
        }

        #read in contents and populate array
        open(my $fileHandle, "<", $file) or die "Can't open file $file: $!\n";
        # skip header
        readline $fileHandle;

        my $index = 0;
        while (my $line = <$fileHandle>) {
            chomp $line;
            $arrayRef->[$index] += $line;
            $index += 1;
        }

        close $fileHandle;
        #print "$prophage\n";
        #print scalar @{$arrayRef} . "\n\n";
    }
    else {
        print "\nEntering dir: $file\n" if $verbose;
    }
}

#Iterate through all prophages encountered, printing their combined index 
#"charts" to separate files in the output directory
sub print_combined_charts {
    my @hashKeys = sort keys %phageHash;
    foreach my $key (@hashKeys) {
        my $outputPath = "$outputDir/$key" . "Chart.txt";
        checkFile($outputPath);
        open(my $outputHandle, ">", $outputPath) or die "Can't open $outputPath: $!";

        my @indexArray = @{$phageHash{$key}};
        my $contents = "$key\n";

        foreach my $index (@indexArray) {
            $contents .= "$index\n";
        }

        #print "$key\n";
        #print scalar @indexArray . "\n\n";

        print $outputHandle "$contents";
        close $outputHandle;
    }

}

sub checkFile {
    my $file = $_[0];

    if (-f $file && !$force) {
        die "$file already exists! To overwrite any files that already exist, rerun with --force option.\n";
  }
}

sub help {
    #expects everything to be a file or directory!
    print "
# combine_charts.pl: Scans through a directory of nucleotide index directories, 
# combining their nucleotide counts
# expects all objects encountered in the input directory to be either files or 
# directories

Options:
    --inputdir: Path to input directory. Expects all items in this directory to 
        be either files or directories.
    --outputdir: Path to directory intended for total nucleotide counts.
    --verbose: Prints updates as the program iterates through inputdir.
    --force: If output files already exist, overwrite them.
    --help: Print this help page.
\n";
}
