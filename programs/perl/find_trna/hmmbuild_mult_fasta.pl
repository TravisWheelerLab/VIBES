#!/usr/bin/env perl

# runs hmmbuild on all .afa files in a directory, outputing resultant hmms in output dir 
use strict;
use warnings;

use Getopt::Long;

# Command-line options
my $inputDir = '';
my $outputPath = '';
my $prefix = '';
my $cpu = 0;
my $help = '';
my $verbose = 0;
my $force = '';

GetOptions (
    "input_dir=s"   => \$inputDir,
    "output=s"      => \$outputPath,
    "prefix=s"      => \$prefix,
    "cpu=i"         => \$cpu,
    "help"          => \$help,
    "verbose"       => \$verbose,
    "force"         => \$force
    )
or die("Unknown option, try --help\n");

if ($help) {
    help();
}
else {
    die "You haven't finished fixing me yet, dummy\n";
    # store all files in input directory in array. We expect that all files in input dir are .afa alignments
    my @afas = glob "$inputDir/*";
    my @outputs = [];

    foreach my $afa (@afas) {
        # match as many characters as possible until we encounter the last '/', then grab everything until we come across a '.' character.
        # This should capture the name of the .afa from its file path
        $afa =~ /.+\/(.+?)\./;
        my $afaName = $1;

        # if prefix has been specified, stick an '_' character on the end
        if ($prefix) {
            $prefix = "$prefix\_";
        }

        # run hmmbuild, then hmmpress
        
        print $outputHandle do_cmd("hmmbuild --cpu $cpu $outputPath") . "\n";
    }
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

sub help {
    print("
# $0: Automatically runs hmmbuild and hmmpress on all .afas in a directory. 
#The program uses the --auto flag to automatically set MAFFT options.
    Usage: perl $0 --input_dir /path/ --output_dir /path/

    Input:
        --input_dir /path/: Path to directory of input .afa files.

    Output:
        --output_dir /path/: Path to directory where output hmms will be stored.
        Hmmpress is also run by $0, so related indexing files will be stored in
        the output dir.
        --prefix str: Prefix that will be appended to output file names.

    Misc:
        --cpu i: number of parallel workers to use in hmmbuild
        --help: Prints this help page.
        --force: Overwrite output files, if they already exist.
        --verbose: Prints commands run by script and their results.\n");
}