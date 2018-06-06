#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use File::Temp qw(tempfile);
use Cwd;
#use Dfamscan;

my $VERSION = 1.4;

$|++;

my $logger_missing;

BEGIN {
  # we need to check if Log4perl is installed
  # if not then we cant use logdie;
  eval {
    require Log::Log4perl;
    Log::Log4perl->import(':easy');
  };

  if ($@) {
    $logger_missing = 1;
  }
};


my $e_max = 10000;
my $score_dominate_default  = 5;
my $rph_trim_default        = 10;
my $min_cov_frac_default    = 0.75;
my $default_cpu             = 8;

main();

1;

sub main {
   my $args = processCommandLineArgs();

   verify_programs($args);
   my ($hits, $header) = get_nhmmscan_hits($args); #array pointer
   my $cnt_before = 1 + $#{$hits};


   if ($args->{no_rph_removal} ) {
      printf STDERR ("nhmmer hits (no overlap removal):   %7d\n", $cnt_before);
   } else {
      $hits = Dfamscan::filter_covered_hits($hits, $args);
      logdie("Error running covered-hit-filter\n") unless ($hits);
      my $cnt_after = 1 + $#{$hits};
      printf STDERR ("nhmmer hits before overlap removal: %7d\n", $cnt_before);
      printf STDERR ("nhmmer hits after overlap removal:  %7d\n", $cnt_after);
   }


   my $sorted = [];
   if ( $args->{sort_method} eq "seq") {
      $sorted = Dfamscan::by_seqandorientandpos($hits);
   } elsif ($args->{sort_method} eq "model") {
      $sorted = Dfamscan::by_modelandseqandorientandpos($hits);
   } elsif ($args->{sort_method} eq "eval") {
      $sorted = Dfamscan::by_evalue($hits);
   }

   open DFOUT, ">$args->{dfam_outfile}" or logdie ("Can't open $args->{dfam_outfile}: $!");
   print DFOUT $header;
   foreach my $row (@$sorted) {
     print DFOUT $row->{line};
   }
   close DFOUT;

   if ($args->{trf_outfile}) {
      my $cnt = run_trf($args->{fastafile}, $args->{trf_outfile});
      logdie("Error running TRF\n") if ($cnt == -1);
      printf STDERR ("TRF repeat count:                   %7d\n", $cnt);
   }

}

sub logdie {
  my $message = shift;
  if ($logger_missing) {
    die $message;
  }
  else {
    LOGDIE($message);
  }
}

sub processCommandLineArgs {

  my %args = ( cpu => $default_cpu, );

  &GetOptions( \%args,
    "help",
    "version",
    "dfam_infile=s",
    "fastafile=s",
    "hmmfile=s",
    "dfam_outfile=s",
    "E=f",
    "T=f",
    "masking_thresh|cut_ga",
    "annotation_thresh|cut_tc",
    "species=i",
    "trf_outfile=s",
    "cpu=i",
    "no_rph_removal",
    "rph_trim=i",
    "score_dominate=f",
    "min_cov_frac=f",
    "sortby_model",
    "sortby_seq",
    "sortby_eval",
    "log_file=s",
  )
  or logdie ("Unknown option, try running --help for more infortmation.");

  help() if ($args{help});
  version() if ($args{version});

  if ($args{log_file}) {
    if ($logger_missing) {
      warn "Log::Log4perl failed to load. No logs will be created\n";
    }
    else {
      if ( $args{log_file} =~ m[^(.+)/]) { # some other directory, it better exist
        die("The directory $1 does not exist")  unless (-d $1);
      }
      Log::Log4perl->easy_init( {  file    => ">$args{log_file}" } );
    }
  }

  if ($args{dfam_infile}) { #queries have already been run
    help("Illegal flags used in addition to --dfam_infile" ) if ($args{fastafile} || $args{hmmfile} || $args{trf_outfile} || $args{masking_thresh} || $args{annotation_thresh} || $args{species});
    logdie("Unable to open $args{dfam_infile}: $!")  unless -e $args{dfam_infile};
  }
  elsif ( $args{fastafile} && $args{hmmfile} ) { #need to run nhmmer
    help("Flag --dfam_infile may not be used with --hmmfile") if ($args{dfam_infile});
    logdie("Unable to open $args{fastafile}: $!")  unless -e $args{fastafile};
    logdie("Unable to open $args{hmmfile}: $!")    unless -e $args{hmmfile};
  }
  else {
    help("Use either (--dfam_infile) or (--fastafile and --hmmfile)");
  }

  if ( $args{dfam_outfile} ) {
    # does the containing directory exist?
    if ( $args{dfam_outfile} =~ m[^(.+)/]) { #some other directory, it better exist
      logdie("The directory $1 does not exist")  unless (-d $1);
    }
   }
   else {
     help("Must specify --dfam_outfile");
   }

   if ( $args{trf_outfile} ) {
     # does the containing directory exist?
     if ( $args{trf_outfile} =~ m[^(.+)/]) { #some other directory, it better exist
       logdie("The directory $1 does not exist")  unless (-d $1);
     }
   }


   my $thresh_arg_cnt=0;
   $args{mask_method} = "TC"; #default
   if ($args{E}) {
      $thresh_arg_cnt++;
      $args{mask_method} = "E";
      logdie("Invalid E value: $!")    if ( $args{E} <= 0 || $args{E} > $e_max ) ;
   } elsif ($args{T}) {
      $thresh_arg_cnt++;
      $args{mask_method} = "T";
   } elsif ($args{masking_thresh}) {
      $thresh_arg_cnt++;
      $args{mask_method} = "GA";
   } elsif ($args{annotation_thresh}) {
      $thresh_arg_cnt++;
      $args{mask_method} = "TC";
   } elsif ($args{species}) {
      $thresh_arg_cnt++;
      $args{mask_method} = "SP";
   }
   if ($thresh_arg_cnt > 1) { #only one should be specified
      help("Only one threshold method should be given");
   }


   my $sort_arg_cnt=0;
   $args{sort_method} = "seq"; #default
   if ($args{sortby_eval}) {
      $sort_arg_cnt++;
      $args{sort_method} = "eval";
   } elsif ($args{sortby_model}) {
      $sort_arg_cnt++;
      $args{sort_method} = "model";
   } elsif ($args{sortby_seq}) {
      $sort_arg_cnt++;
      $args{sort_method} = "seq";
   }
   if ($sort_arg_cnt > 1) { #only one should be specified
      help("Only one sort method should be specified");
   }

   if ($args{no_rph_removal}) {
     if ($args{rph_trim}) {
        help("Don't specify both --no_rph_removal and --rph_trim" );
     }
   } else {
      $args{rph_trim} = $rph_trim_default;
   }


   unless ($args{score_dominate}) {
      $args{score_dominate} = $score_dominate_default;
   }

   unless ($args{min_cov_frac}) {
      $args{min_cov_frac} = $min_cov_frac_default;
   }


   return(\%args);
}


sub version {
   my $args = shift;
   printf ("%16s : version %s\n", $0 , $VERSION);

   if (in_path('nhmmscan')) {
      my $res = `nhmmscan -h`;
      $res =~ /^# (HMMER .+?);/m;
      printf ("%16s : version %s\n", "nhmmscan" , $1);
   }

   if (in_path('trf')) {
      my $res = `trf`;
      $res =~ /Tandem Repeats Finder, Version (\S+)/m;
      printf ("%16s : version %s\n", "trf", $1);
   }
   exit;
}

sub help {

print STDERR <<EOF;

Command line options for controlling $0
-------------------------------------------------------------------------------

   --help       : prints this help messeage
   --version    : prints version information for this program and
                  both nhmmscan and trf


   Requires either
    --dfam_infile <s>    Use this is you've already run nhmmscan, and
                         just want to perfom dfamscan filtering/sorting.
                         The file must be the one produced by nhmmscan's
                         --dfamtblout flag.
                         (Note: must be nhmmscan output, not nhmmer output)
   or both of these
    --fastafile <s>      Use these if you want dfamscan to control a
    --hmmfile <s>        run of nhmmscan, then do filtering/sorting

   Requires
    --dfam_outfile <s>   Output file, also in dfamtblout format

   Optionally, one of these  (only -E and -T allowed with --dfam_infile)
    -E <f>               >0, <=$e_max
    -T <f>
    --masking_thresh/--cut_ga
    --annotation_thresh/--cut_tc  Default
    --species <i>        (not yet implemented)

   Optionally one of these
    --sortby_eval
    --sortby_model
    --sortby_seq         Default

   Redundant Profile Hit (RPH) removal (only if not --no_rph_removal)
    --rph_trim <i>       If hit A has higher score than hit B and covers
                         all but at most rph_trim bases on either side,
                         then A dominates B. Default $rph_trim_default
    --score_dominate <f> If hit A has score this much higher (bits) than
                         hit B, then A can dominate B even if it's is
                         shorter, so long as the --min_cov_frac threshold
                         is met. Default $score_dominate_default
    --min_cov_frac <f>   If hit A exceeds hit B by --score_dominate bits,
                         and either A or B is covered at least --min_cov_frac
                         by the other, then A dominates B. Default $min_cov_frac_default

   All optional
    --trf_outfile <s>    Runs trf, put results in <s>; only with --fastafile
    --cpu <i>            Default $default_cpu
    --no_rph_removal     Don't remove redundant profile hits
    --log_file <s>
EOF


  logdie ("\n**********\n\n$_[0]\n") if $_[0];


  exit(1);
}


sub verify_programs {
   my ($args) = @_;

   logdie("nhmmscan not found in \$PATH")  unless (in_path('nhmmscan'));

   if ($args->{trf_outfile}) {
      logdie("trf not found in \$PATH")  unless (in_path('trf'));
   }
}

sub in_path {
   my $prog = $_[0];
   my @dirs = split(/:/,$ENV{'PATH'});
   foreach my $dir (@dirs) {
      #from perl cookbook, 2nd edition
      $dir =~   s{ ^ ~ ( [^/]* ) }
                     { $1
                           ? (getpwnam($1))[7]
                           : ( $ENV{HOME} || $ENV{LOGDIR}
                                || (getpwuid($<))[7]
                             )
                     }ex;

      if (-e "$dir/$prog")  {
         return 1;
      }
   }
   return 0;
}

sub run_trf {
   my ($fastafile, $trf_outfile) = @_;

   # change dir to the input directory as trf wont let us
   # specify where the output should go
   my $cur_dir = getcwd;
   my $out_dir = $fastafile;
   $out_dir =~ s|[^/]*$||;
   chdir($out_dir);
   # params to search for old simple repeats:  2 3 5 75 20 33 7
   my $cmd = "trf $fastafile 2 7 7 80 10 50 10 -d -h > /dev/null";
   my $res = system($cmd);
   rename "$fastafile.2.7.7.80.10.50.10.dat", $trf_outfile or return -1;

   #get counts
   open FH, "<$trf_outfile";
   my $cnt = 0;
   while (my $line = <FH>) {
      $cnt++ if ($line =~ /^\d/);
   }
   close FH;
   chdir($cur_dir);
   return $cnt;
}

sub get_nhmmscan_hits {
   my ($args) = @_;

   my @hits;
   my $header= "";
   my $header_done = 0;
   if ($args->{dfam_infile}) {
      open FH, "<$args->{dfam_infile}";
      while (my $line = <FH>) {
         if ($line =~ /^#/) {
            $header .= $line if  !$header_done;

            $header_done = 1 if ($line =~ /-------------------/); # that'll be the final line of the first header, no need to replicate more headers
            next;
         }
         $header_done = 1;
         if ($args->{mask_method} =~ /^(E|T)$/ ) { #filter hits
            my @vals = split(/\s+/, $line);
            if ($args->{mask_method} eq "E") {
               next if ($vals[4] > $args->{E});
            } elsif ($args->{mask_method} eq "T") {
               next if ($vals[3] < $args->{T});
            }
         }
         push @hits, get_hit_from_hitline($line);

      }
      close FH;
   } else { # ($args{fastafile} && $args{hmmfile});
      my $cmd = "nhmmscan --noali";

      if ($args->{mask_method} eq "SP") {
         #um ...   do I have to split the hmm file in half, part for matching species, and part for not?
         logdie("Don't know how to handle --species yet");
         # two nhmmer jobs?

      } else {

         if ($args->{mask_method} eq "E") {
            $cmd .= " -E $args->{E}";
         } elsif ($args->{mask_method} eq "T") {
            $cmd .= " -T $args->{T}";
         } elsif ($args->{mask_method} eq "GA") {
            $cmd .= " --cut_ga";
         } elsif ($args->{mask_method} eq "TC") {
            $cmd .= " --cut_tc";
         }


         my ($tmpfh, $dfout_filename) = tempfile();

         $cmd .= " --dfamtblout $dfout_filename";
         $cmd .= " --cpu=$args->{cpu}";
         $cmd .= " $args->{hmmfile}";
         $cmd .= " $args->{fastafile}";

         #print ("$cmd\n");
         my $result = system ("$cmd > /dev/null");
         logdie("Error running command:\n$cmd\n") if $result;
         open FH, "<$dfout_filename";
         while (my $line = <FH>) {
            if ($line =~ /^#/) {
               $header .= $line if  !$header_done;
               next;
            }
            $header_done = 1;
            push @hits, get_hit_from_hitline($line);
         }
         close FH;
         unlink $dfout_filename;
      }
   }
   return \@hits, $header;
}






sub get_hit_from_hitline {
   my ($model, $acc, $seq, $score, $eval, $tmp1, $tmp2, $tmp3, $orient, $start, $end) = split(/\s+/,$_[0]);
   if ($orient eq "-") {
      $tmp1 = $start;
      $start = $end;
      $end = $tmp1;
   }

   return {
      model  => $model,
      acc    => $acc,
      seq    => $seq,
      score  => $score,
      e_val  => $eval,
      orient => $orient,
      start  => $start,
      end    => $end,
      line   => $_[0],
    };
}

package Dfamscan;

use strict;

#our $VERSION = 1.4;

$|++;



sub filter_covered_hits {
   my @sorted = @{ by_seqandpos($_[0])};

   my $args   = $_[1];
   my $i=0;
   my $j;
   my $tmp;

   my @dominated;
   my $cur_seq = "";

   my ($model_a, $acc_a, $seq_a, $score_a, $eval_a, $orient_a, $start_a, $end_a);
   my ($model_b, $acc_b, $seq_b, $score_b, $eval_b, $orient_b, $start_b, $end_b);

   for ($i=0; $i<$#sorted; $i++) {
      ($model_a, $acc_a, $seq_a, $score_a, $eval_a, $orient_a, $start_a, $end_a) = get_vals_from_hit($sorted[$i]);

      if ($end_a <= $start_a + 5) {
         $dominated[$i] = 1;
         next;
      }
      #set these to jumpstart the loop below
      $j=$i;
      $seq_b   = $seq_a;
      $start_b = $start_a;
      $end_b   = $end_a;

      if ($seq_a ne $cur_seq) {
         $cur_seq = $seq_a;
      }

      #grab all the hits that are even just-a-little overlapping
      while ($j<$#sorted && $seq_b eq $seq_a && $start_b < $end_a - 5) {
         $end_a = ($end_b > $end_a ? $end_b : $end_a);
         $j++;
         ($model_b, $acc_b, $seq_b, $score_b, $eval_b, $orient_b, $start_b, $end_b) = get_vals_from_hit($sorted[$j]);
      }
      if ($j == $#sorted && $seq_b eq $seq_a && $start_b < $end_a) {
         mask_dominated($i, $j, \@sorted, \@dominated, $args);
         $i = $j; # so the next loop, the first entry will be the one that comes right after the end of the current range.
      } else {
         mask_dominated($i, $j-1, \@sorted, \@dominated, $args);
         $i = $j-1; # so the next loop, the first entry will be the one that comes right after the end of the current range.
      }
   }



   $i = 0;
   $j = 0;
   while ($j <= $#sorted) {
      unless ($dominated[$j]) {
         $sorted[$i] = $sorted[$j];
         $i++;
      }
      $j++;
   }

   #resize the array, to drop off leftover cruft;
   $#sorted = $i-1;

   return \@sorted;

}

sub mask_dominated {
   my ($first, $last, $sorted_ref, $dominated_ref, $args) = @_;
   my $i;
   my $j;

   return if $first == $last; # no chance one hit dominates itself

   #this needs to get the list of un-dominated hits, then find the ones that
   #those dominate, then repeat with the remainder
   my $todo_cnt = $last-$first+1;

   my @todo_mask;
   for ($i=$first; $i<=$last; $i++) {
      $todo_mask[$i] = 1;
   }

   my @unbeaten;
   while ($todo_cnt > 0) {
      my %beaten;
      for ($i=$first; $i<=$last-1; $i++) {
         next if ( !$todo_mask[$i] || $beaten{$i});

         # while going through the list to find hits that are not dominated,
         # keep a hash of, for each hit, the ones it dominates
         for ($j=$i+1; $j<=$last; $j++) {
            next if ( !$todo_mask[$j] || $beaten{$j});
            my $res = compare_hits ( $$sorted_ref[$i], $$sorted_ref[$j], $args);
            if ($res == 1) {
               $beaten{$j} = 1;
            } elsif ($res == -1) {
               $beaten{$i} = 1;
               last; # i has been dominated, so stop
            }
         }
      }

      for ($i=$first; $i<=$last; $i++) {
         next if ( !$todo_mask[$i] || $beaten{$i});
         for ($j=$first; $j<=$last; $j++) {
            next if $i == $j;
            next if !$todo_mask[$j];
            # Since we didn't track all pairs of hits, just re-test which hits are dominated by an unbeaten hit.
            my $res = compare_hits ( $$sorted_ref[$i], $$sorted_ref[$j], $args);

            if ($res == 1) {
               $todo_cnt--;
               $todo_mask[$j] = 0;
               $$dominated_ref[$j] = 1;
            }
         }
         push (@unbeaten, $i);
         $todo_mask[$i] = 0;
         $todo_cnt--;
      }
   }

   @unbeaten = sort {$a<=>$b} @unbeaten;

}


#if A dominates, return 1;  if B dominates, return -1, else return 0
sub compare_hits {
   my $tmp;
   my ($hit_a, $hit_b, $args) = @_;
   my ($i, $j);

   my ($covered_start, $covered_end, $covered_len, $a_len, $b_len);


   my $swap_mod=1;
   if ($hit_b->{start} < $hit_a->{start} || ($hit_b->{start} == $hit_a->{start} && $hit_b->{end} > $hit_a->{end}) ) {
      my $tmp = $hit_a;
      $hit_a = $hit_b;
      $hit_b = $tmp;
      $swap_mod = -1;
   }


   if ( $hit_a->{score} == $hit_b->{score} ) {
      if ( $hit_a->{end}-$hit_a->{start} >= $hit_b->{end}-$hit_b->{start}) {
          return $swap_mod;
      } else {
          return -1 * $swap_mod;
      }
   } elsif ( $hit_a->{score} >= $hit_b->{score} &&  $hit_a->{end} >= $hit_b->{end} - $args->{rph_trim}) { # A is clearly superior both in score and (essentially) in boundary
      return $swap_mod;
   } elsif ( $hit_a->{score} >= $hit_b->{score} + $args->{score_dominate} ) {
      # A has (default ~100x) better E-value; if one of the two hits is mostly covered by the other, A dominates B
      $covered_start = $hit_a->{start} > $hit_b->{start} ? $hit_a->{start} : $hit_b->{start};
      $covered_end   = $hit_a->{end} < $hit_b->{end}     ? $hit_a->{end}   : $hit_b->{end};
      $covered_len   = $covered_end-$covered_start +1;
      $a_len = $hit_a->{end} - $hit_a->{start} + 1;
      $b_len = $hit_b->{end} - $hit_b->{start} + 1;
      if ( $covered_len/$a_len >= $args->{min_cov_frac} || $covered_len/$b_len >= $args->{min_cov_frac} ) {
         return $swap_mod;
      }
   } elsif ( $hit_b->{score} >= $hit_a->{score} &&  $hit_b->{start} <= $hit_a->{start} + $args->{rph_trim} ) { # B is cleanly superior both in score and (essentially) in boundary
      return -1 * $swap_mod;
   } elsif ( $hit_b->{score} >= $hit_a->{score} + $args->{score_dominate}) {
      # B has (default ~100x) better E-value; if one of the two hits is mostly covered by the other, B dominates A
      $covered_start = $hit_a->{start} > $hit_b->{start} ? $hit_a->{start} : $hit_b->{start};
      $covered_end   = $hit_a->{end} < $hit_b->{end}     ? $hit_a->{end}   : $hit_b->{end};
      $covered_len   = $covered_end-$covered_start +1;
      $a_len = $hit_a->{end} - $hit_a->{start} + 1;
      $b_len = $hit_b->{end} - $hit_b->{start} + 1;
	  if ( $covered_len/$a_len >= $args->{min_cov_frac} || $covered_len/$b_len >= $args->{min_cov_frac} ) {
         return -1 * $swap_mod;
      }
   }

   return 0;

}


sub get_vals_from_hit {
   my $a      = $_[0];

   my $model  = $a->{model};
   my $acc    = $a->{acc};
   my $seq    = $a->{seq};
   my $score  = $a->{score};
   my $eval   = $a->{e_val};
   my $orient = $a->{orient};
   my $start  = $a->{start};
   my $end    = $a->{end};

   return ($model, $acc, $seq, $score, $eval, $orient, $start, $end);

}


## sorting helper functions

sub by_modelandseqandorientandpos {
  my $input = shift;
  my @result = sort {
    ($a->{model} cmp $b->{model})
    || ($a->{seq} cmp $b->{seq})
    || ($a->{orient} cmp $b->{orient})
    || ($a->{start} <=> $b->{start})
    || ($a->{end}   <=> $b->{end})
    || ($a->{score} <=> $b->{score})
  } @$input;
  return \@result;
}


sub by_seqandpos {
  my $input = shift;
  my @result = sort {
    ($a->{seq} cmp $b->{seq})
    || ($a->{start} <=> $b->{start})
    || ($a->{end}   <=> $b->{end})
    || ($a->{score} <=> $b->{score})
    || ($a->{model} cmp $b->{model})
  } @$input;
  return \@result;
}

sub by_seqandorientandpos {
  my $input = shift;
  my @result = sort {
    ($a->{seq} cmp $b->{seq})
    || ($a->{orient} cmp $b->{orient})
    || ($a->{start} <=> $b->{start})
    || ($a->{end}   <=> $b->{end})
    || ($a->{score} <=> $b->{score})
    || ($a->{model} cmp $b->{model})
  } @$input;
  return \@result;
}

sub by_model {
  my $input = shift;
  my @result = sort {
    ($a->{model} cmp $b->{model})
    || ($a->{seq} cmp $b->{seq})
    || ($a->{orient} cmp $b->{orient})
    || ($a->{start} <=> $b->{start})
    || ($a->{score} <=> $b->{score})
  } @$input;
  return \@result;
}

sub by_evalue {
  my $input = shift;
  my @result = sort {
    ($a->{e_val} <=> $b->{e_val})
    || ($a->{score} <=> $b->{score})
    || ($a->{seq} cmp $b->{seq})
    || ($a->{orient} cmp $b->{orient})
    || ($a->{start} <=> $b->{start})
    || ($a->{model} cmp $b->{model})
  } @$input;
  return \@result;
}

sub by_score {
  my $input = shift;
  my @result = sort {
    ($a->{score} <=> $b->{score})
    || ($a->{e_val} <=> $b->{e_val})
    || ($a->{seq} cmp $b->{seq})
    || ($a->{orient} cmp $b->{orient})
    || ($a->{start} <=> $b->{start})
    || ($a->{model} cmp $b->{model})
  } @$input;
  return \@result;
}

1;
