#!/usr/bin/env perl

# runs mafft on all .fasta files in a directory, outputing resultant aligned fastas in output dir 
use strict;
use warnings;

use Getopt::Long;

# Command-line options
my $inputDir = '';
my $outputDir = '';
my $prefix = '';
my $thread = 0;
my $help = '';
my $verbose = 0;
my $force = '';

GetOptions (
    "input_dir=s"   => \$inputDir,
    "output_dir=s"  => \$outputDir,
    "prefix=s"      => \$prefix,
    "thread=i"      => \$thread,
    "help"          => \$help,
    "verbose"       => \$verbose,
    "force"         => \$force
    )
or die("Unknown option, try --help\n");

if ($help) {
    help();
}
else {
    # store all files in input directory in array. We expect that all files in input dir are genomes
    my @fastas = glob "$inputDir/*";

    foreach my $fasta (@fastas) {

        # open .fasta file and verify it has at least 2 entries
        my $entryCount = 0;
        open(my $fastaHandle, "<", $fasta) or die "Can't open fasta file $fasta: $!\n";
        while (my $line = <$fastaHandle>) {
            if ($line =~ /^\>/) {
                $entryCount += 1;

                # if we have at least 2 entries, we can safely run MAFFT, so we exit loop
                if ($entryCount > 1) {
                    last;
                }
            }
        }

        close($fastaHandle);

        # if at least 2 fasta entries are in the .fasta file, we can safely run MAFFT
        if ($entryCount > 1) {
            # match as many characters as possible until we encounter the last '/', then grab everything until we come across a '.' character.
            # This should capture the name of the genome from its path
            $fasta =~ /.+\/(.+?)\./;
            my $fastaName = $1;

            # if prefix has been specified, stick an '_' character on the end
            if ($prefix) {
                $prefix = "$prefix\_";
            }

            # run MAFFT with --maxiterate 1000 to increase sensitivity 
            do_cmd("mafft --thread $thread --quiet --maxiterate 1000 $fasta > $outputDir/$prefix$fastaName.afa");
        }
        else {
            print STDERR "Warning: $fasta contains only 1 entry. No .afa file will be generated!\n";
        }
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
# $0: Automatically runs MAFFT on all .fastas in a directory. 
#The program uses the --auto flag to automatically set MAFFT options.
    Usage: perl $0 --input_dir /path/ --output_dir /path/

    Input:
        --input_dir /path/: Path to directory of target genomes in .fasta
            format. All files in this directory are expected to be genomes.
            This option is mutually exclusive with --input_fasta.

    Output:
        --output_dir /path/: Path to directory where output tRNA sequences
            will be stored in .fasta format. Since these files are intended
            to generate a MSA, tRNA sequences are saved in files named after
            the respective amino acids they encode for.
        --prefix str: Prefix that will be appended to output file names.

    Misc:
        --help: Prints this help page.
        --force: Overwrite output files, if they already exist.
        --verbose: Prints commands run by script and their results.\n");
}