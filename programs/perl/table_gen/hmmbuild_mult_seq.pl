#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $inputPath = '';
my $outputPath = '';
my $amino = 0;
my $dna = 0;
my $rna = 0;
my $verbose = 0;
my $help = 0;

GetOptions (
    "input=s"   => \$inputPath,
    "output=s"  => \$outputPath,
    "amino"     => \$amino,
    "dna"       => \$dna,
    "rna"       => \$rna,
    "verbose"   => \$verbose,
    "help"      => \$help
    )
or die("Unknown argument, try --help\n");

if($help) {
    help();
    exit
}

if (($dna + $rna + $amino) > 1) {
    print "Warning: only one of --dna, --rna, or --amino can be used at the same time\n";
    exit;
}

my $pathToFolder = $inputPath;
$pathToFolder =~ s/[\w-]+.fasta//; # remove file name from path to .fasta file, giving you the path to its folder

# open .fasta file, demarcate lines with > rather than /n
open(my $fileHandle, "<", $inputPath) or die "Can't open .fasta file $: $!";
{
    local $/ = "\n>"; # set end of input record separator to \n>, so it splits on lines starting with >. Unfortunately,
    # this is necessary because very rarely, .fasta header lines will contain > characters somewhere in the middle
    my $inc = 1;

    # grab each entry, save it in a temporary file, run it through hmmbuild, and delete it
    while (my $entry = <$fileHandle>)  {

        # remove any > characters at beginning or end of entry, leaving any in the middle of header lines intact
        $entry =~ s/>$//;
        $entry =~ s/^>//;
        # capture header, sequence data separately
        $entry =~ m/(.+?)\n(.+)/s;
        my $name = $1;
        my $seq = $2;
        # remove any whitespace in sequence data
        $seq =~ s/\s//g;


        my $tempFastaFile = "temp$inc";
        
        $tempFastaFile = "$pathToFolder/$tempFastaFile.fasta"; #create file path to temporary .fasta file


        open(my $fastaFile, ">", $tempFastaFile) or die "Can't create temporary .fasta file at $tempFastaFile: $!";
        {
            print $fastaFile ">$name\n$seq";
        }

        close $fastaFile;

        my $tempHmmFile = $tempFastaFile;
        $tempHmmFile =~ s/\.fasta/.hmm/; # replace .fasta with .hmm in file path, to create temporary .hmm file with same name

        # create command, depending on options specified at command line
        my $hmmbuildCmd = "hmmbuild ";

        if($dna) {
            $hmmbuildCmd .= "--dna ";
        }
        elsif($rna) {
            $hmmbuildCmd .= "--rna ";
        }
        elsif($amino) {
            $hmmbuildCmd .= "--amino ";
        }

        $hmmbuildCmd .= "$tempHmmFile $tempFastaFile";

        do_cmd($hmmbuildCmd);

        if($verbose) {
            print "hmmbuild round $inc completed\n\n";
        }

        $inc++;
    }
}

close $fileHandle;

do_cmd("cat $pathToFolder/temp*.hmm > $outputPath");
do_cmd("rm $pathToFolder/temp*");

sub do_cmd {
    my $cmd = $_[0];

    if ($verbose > 0) {
        print "$cmd\n";
    }

    return `$cmd`;
}