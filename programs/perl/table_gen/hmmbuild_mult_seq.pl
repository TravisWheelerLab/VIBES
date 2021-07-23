#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $inputPath = '';
my $outputPath = '';
my $cpuCount = 0;
my $amino = 0;
my $dna = 0;
my $rna = 0;
my $verbose = 0;
my $help = 0;

GetOptions (
    "input=s"   => \$inputPath,
    "output=s"  => \$outputPath,
    "cpu=i"     => \$cpuCount,
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

unless (-f $inputPath){
    die "Input $inputPath is not a plain file: $!"
}

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
        $entry =~ s/^([^\w]+?>|>)//;
        # capture header, sequence data separately
        $entry =~ m/(.+?)\n(.+)/s;
        # we add a newline to the header to ensure it ends with a whitespace character, which makes a line of Regex complain less
        my $header = "$1\n";
        my $seq = $2;
        # remove any non-word characters in sequence data
        $seq =~ s/[^\w]//g;

        my $tempFastaFile = "$inputPath.temp.$inc";

        open(my $fastaFile, ">", $tempFastaFile) or die "Can't create temporary .fasta file at $tempFastaFile: $!";
        {   
            # no additional newline required- we added one above
            print $fastaFile ">$header$seq";
        }

        close $fastaFile;

        my $tempHmmFile = "$inputPath.hmm.temp.$inc";

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

        # set number of cpus for hmmbuild to use
        if ($cpuCount > 0) {
            $hmmbuildCmd .= "--cpu $cpuCount ";
        }

        my $name = "";

        # grab header line up to first whitespace character
        if ($header =~ m/(.+?)\s/) {
            $name = $1;
        }
        else {
            die "Can't parse $header\n";
        }

        # set name of sequence to header line
        $hmmbuildCmd .= "-n \"$name\" ";

        $hmmbuildCmd .= "$tempHmmFile $tempFastaFile";

        # temporary bugfixing print
       # my $path = $ENV{'PATH'};
       # print "PATH: $path\n";
       # my $result = `type hmmbuild`;
       # print "HMMBUILD: $result\n";

        do_cmd($hmmbuildCmd);

        if($verbose) {
            print "hmmbuild round $inc completed\n\n";
        }

        $inc++;
    }
}

sub help {
    print "
#$0
#Runs hmmbuild on each entry in input .fasta file, creating a HMM for each.
#Automatically concatenates all HMMS into one output file
-------------------------------------------------------------------------------
--help: Prints this message

Required:
 --input <s>        Input .fasta file of multiple sequences we want to create
                    HMMS for
 --output <s>       Output .hmm file to store resultant HMMS in

Optional:
 --verbose          Prints information about commands used, how many .fasta
                    entries have been hmmbuilt
 --cpu <i>          How many threads hmmbuild will use (i > 0)

 Only one of the following 3 optional flags can be used at a time:
  --amino           Specifies that .fasta entries contain amino acid seq
  --dna             Specifies that .fasta entries contain DNA seq
  --rna             Specifies that .fasta entries contain RNA seq
 
";
}

close $fileHandle;

do_cmd("cat $inputPath.hmm.temp.* > $outputPath");
do_cmd("rm $inputPath.temp.*");
do_cmd("rm $inputPath.hmm.temp.*");

do_cmd("hmmpress $outputPath");

sub do_cmd {
    my $cmd = $_[0];

    if ($verbose > 0) {
        print "$cmd\n";
    }

    return `$cmd`;
}