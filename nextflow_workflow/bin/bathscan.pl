#!/usr/bin/env perl

#derived from dfamscan.pl

use strict;
use warnings;
use Getopt::Long;
use File::Temp qw(tempfile);
use Cwd;
use File::Basename;

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
my $min_cov_frac_default    = 0.75;

main();

1;

sub main {
   my $args = processCommandLineArgs();


   my ($hits, $header) = get_frahmmscan_hits($args); #array pointer
   my $cnt_before = 1 + $#{$hits};


   if ($args->{no_rph_removal} ) {
      printf STDERR ("nhmmer hits (no overlap removal):   %7d\n", $cnt_before);
   } else {
      $hits = &filter_covered_hits_using_masklevel($hits, $args->{min_cov_frac}, 1);
      logdie("Error running covered-hit-filter\n") unless ($hits);
      my $cnt_after = 1 + $#{$hits};
      printf STDERR ("nhmmer hits before overlap removal: %7d\n", $cnt_before);
      printf STDERR ("nhmmer hits after overlap removal:  %7d\n", $cnt_after);
   }


   my $sorted = [];
   if ( $args->{sort_method} eq "seq") {
      $sorted = by_seqandorientandpos_local($hits);
   } elsif ($args->{sort_method} eq "model") {
      $sorted = by_modelandseqandorientandpos_local($hits);
   } elsif ($args->{sort_method} eq "eval") {
      $sorted = by_evalue_local($hits);
   }

   open DFOUT, ">$args->{outfile}" or logdie ("Can't open $args->{outfile}: $!");
   print DFOUT $header;
   foreach my $row (@$sorted) {
     print DFOUT $row->{line};
   }
   close DFOUT;

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

  my %args = ( );

  &GetOptions( \%args,
    "help",
    "version",
    "infile=s",
    "outfile=s",
    "E=f",
    "T=f",
    "masking_thresh|cut_ga",
    "annotation_thresh|cut_tc",
    "no_rph_removal",
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

  logdie("Unable to open $args{infile}: $!")  unless -e $args{infile};

  if ( $args{outfile} ) {
    # does the containing directory exist?
    if ( $args{outfile} =~ m[^(.+)/]) { #some other directory, it better exist
      logdie("The directory $1 does not exist")  unless (-d $1);
    }
   }
   else {
     help("Must specify --outfile");
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


   unless ($args{min_cov_frac}) {
      $args{min_cov_frac} = $min_cov_frac_default;
   }


   return(\%args);
}


sub version {
   my $args = shift;
   printf ("%16s : version %s\n", $0 , $VERSION);


   exit;
}

sub help {

print STDERR <<EOF;
Command line options for controlling $0
-------------------------------------------------------------------------------
   --help       : prints this help messeage
   --version    : prints version information for this program
   Requires
    --infile <s>         Use this to provide the file that holds the
                         frahmmer results that you've already computed
   and
    --outfile <s>   Output file, also in tblout format
   Optionally, one of these  (only -E and -T allowed with --infile)
    -E <f>               >0, <=$e_max
    -T <f>

   Optionally one of these
    --sortby_eval
    --sortby_model
    --sortby_seq         Default
   Redundant Profile Hit (RPH) removal (only if not --no_rph_removal)
    --min_cov_frac <f>   This is similar to the MaskLevel concept in
                         crossmatch.  A match is considered non-redundant
                         if at least (100-min_cov_frac)% of it's aligned
                         bases are not contained within the domain of any
                         higher-scoring hit. Default: $min_cov_frac_default
   All optional
    --no_rph_removal     Don't remove redundant profile hits
    --log_file <s>
EOF


  logdie ("\n**********\n\n$_[0]\n") if $_[0];


  exit(1);
}


sub get_frahmmscan_hits {
   my ($args) = @_;

   my @hits;
   my $header= "";
   my $header_done = 0;
   if (!$args->{infile}) {
      logdie("frahmmscan requires that you first run frahmmer\n")
   }

   open FH, "<$args->{infile}";
   while (my $line = <FH>) {
      if ($line =~ /^#/) {
         $header .= $line if  !$header_done;

         $header_done = 1 if ($line =~ /-------------------/); # that'll be the final line of the first header, no need to replicate more headers
         next;
      }
      $header_done = 1;
      if ($args->{mask_method} =~ /^(E|T)$/ ) { #filter hits
         my @vals = split(/\s+/, $line);
         print "[[ $vals[12] ]]  -> [[ $vals[13] ]]\n";
         if ($args->{mask_method} eq "E") {
            next if ($vals[12] > $args->{E});
         } elsif ($args->{mask_method} eq "T") {
            next if ($vals[13] < $args->{T});
         }
      }
      push @hits, get_hit_from_hitline($line);

   }
   close FH;

   return \@hits, $header;
}


sub get_hit_from_hitline {
    # Target, acc, Query, query acc, hmm len, hmm from, hmm to, ali from, ali to, env from, env to, evalue, score, etc
   my ($seq, $acc, $model, $tmp1, $tmp2, $tmp3, $tmp4, $tmp5, $start, $end, $tmp6, $tmp7, $eval, $score, $tmp8) = split(/\s+/,$_[0]);
   my $orient = "+";
   if ($start > $end) {
      $tmp1 = $start;
      $start = $end;
      $end = $tmp1;
      $orient = "-";
   }

    # TODO: remove me
    # printf("$model\n$acc\n$seq\n$start\n$end\n$eval\n$score\n$orient\n\n\n");

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


sub filter_covered_hits_using_masklevel {
  my $sHits = by_seq_start_end_score_local($_[0]);
  my $masklevel = $_[1];
  my $keepLowScoringSigExt = $_[2];
  my $DEBUG = 0;

  my %deleteHash = ();    # Hash to hold indexes of filtered entries
  my @scores = ();
  my $highest_prev_end = 0;
  my $prev_seq = "";
  for ( my $i = 0 ; $i <= $#{$sHits}; $i++ ) {
    my $hit = $sHits->[$i];
    my $hit_seq = $hit->{seq};
    my $hit_start = $hit->{start};
    my $hit_end = $hit->{end};
    my $hit_score = $hit->{score};

    # If we are still on the same sequence *and* this
    # overlaps something we already have, then simply add it 
    # and continue
    if ( $hit_seq eq $prev_seq && $hit_start < $highest_prev_end - 5 ) {
      push @scores, [ $hit_score, $i, ($hit_end - $hit_start) ];
    }else {  
      # Either we have moved on to another sequence and/or 
      # found a hit that doesn't overlap anything in @scores
      # (or @scores is empty).  Only process clusters with >1
      # annotation ( singletons are not dominated by anything ).
      if ( @scores > 1 ) {
        # Sort by score descending, then by length descending
        @scores = sort { ($b->[0] <=> $a->[0] ) || ($b->[2] <=> $a->[2]) } @scores;
        for ( my $j = 0; $j <= $#scores; $j++ ){
           my $jidx = $scores[$j]->[1];
           if ( exists $deleteHash{$jidx} && $deleteHash{$jidx} == 1 ){
              next;
           }
           for ( my $k = $j+1; $k <= $#scores; $k++ ){
             my $kidx = $scores[$k]->[1];
             if ( exists $deleteHash{$kidx} && $deleteHash{$kidx} == 1 ) {
               next;
             }
             my $hit_a = $sHits->[ $jidx ];
             my $hit_a_start = $hit_a->{start};
             my $hit_a_end = $hit_a->{end};
             my $hit_a_score = $hit_a->{score};
             my $hit_b = $sHits->[ $kidx ];
             my $hit_b_start = $hit_b->{start};
             my $hit_b_end = $hit_b->{end};
             my $hit_b_score = $hit_b->{score};

             # Short circuit if no overlap possible
             if ( $hit_b_end < $hit_a_start || $hit_b_start > $hit_a_end ) {
               next;
             }

             my $covered_start = $hit_a_start;
             $covered_start = $hit_b_start if ($hit_a_start <= $hit_b_start);
             my $covered_end = $hit_a_end;
             $covered_end = $hit_b_end if ($hit_a_end >= $hit_b_end );
             my $covered_len = $covered_end - $covered_start + 1;
             my $a_len = $hit_a_end - $hit_a_start + 1;
             my $b_len = $hit_b_end - $hit_b_start + 1;
             if ( $keepLowScoringSigExt ) {
                if ( ($covered_len / $b_len) >= $masklevel ) {
                    $deleteHash{$kidx} = 1;
                    next;
                }
             }else {
                if ( ($covered_len / $a_len) >= $masklevel ||
                     ($covered_len / $b_len) >= $masklevel ) {
                    $deleteHash{$kidx} = 1;
                    next;
                }
             }
           }
         } # for ( j
       } # if
       @scores = ();
       push @scores, [ $hit_score, $i, ($hit_end - $hit_start) ];
     }
     $prev_seq = $hit_seq;
     if ( $hit_end > $highest_prev_end ) {
         $highest_prev_end = $hit_end;
     }
   }
   #### Trailing case
      if ( @scores > 1 ) {
        # Sort by score descending, then by length descending
        @scores = sort { ($b->[0] <=> $a->[0] ) || ($b->[2] <=> $a->[2]) } @scores;
        for ( my $j = 0; $j <= $#scores; $j++ ){
           my $jidx = $scores[$j]->[1];
           if ( exists $deleteHash{$jidx} && $deleteHash{$jidx} == 1 ){
              next;
           }
           for ( my $k = $j+1; $k <= $#scores; $k++ ){
             my $kidx = $scores[$k]->[1];
             if ( exists $deleteHash{$kidx} && $deleteHash{$kidx} == 1 ) {
               next;
             }
             my $hit_a = $sHits->[ $jidx ];
             my $hit_a_start = $hit_a->{start};
             my $hit_a_end = $hit_a->{end};
             my $hit_a_score = $hit_a->{score};
             my $hit_b = $sHits->[ $kidx ];
             my $hit_b_start = $hit_b->{start};
             my $hit_b_end = $hit_b->{end};
             my $hit_b_score = $hit_b->{score};

             # Short circuit if no overlap possible
             if ( $hit_b_end < $hit_a_start || $hit_b_start > $hit_a_end ) {
               next;
             }

             my $covered_start = $hit_a_start;
             $covered_start = $hit_b_start if ($hit_a_start <= $hit_b_start);
             my $covered_end = $hit_a_end;
             $covered_end = $hit_b_end if ($hit_a_end >= $hit_b_end );
             my $covered_len = $covered_end - $covered_start + 1;
             my $a_len = $hit_a_end - $hit_a_start + 1;
             my $b_len = $hit_b_end - $hit_b_start + 1;
             if ( $keepLowScoringSigExt ) {
                if ( ($covered_len / $b_len) >= $masklevel ) {
                    $deleteHash{$kidx} = 1;
                    next;
                }
             }else {
                if ( ($covered_len / $a_len) >= $masklevel ||
                     ($covered_len / $b_len) >= $masklevel ) {
                    $deleteHash{$kidx} = 1;
                    next;
                }
             }
           }
         } # for ( j
       } # if

  my $i = 0;
  my $j = 0;
  while ($j <= $#{$sHits}) {
     unless ( exists $deleteHash{$j} ) {
        $sHits->[$i] = $sHits->[$j];
        $i++;
     }
     $j++;
  }

  # resize the array, to drop off leftover cruft;
  splice(@{$sHits},$i);

  return $sHits;
}


sub by_seq_start_end_score_local {
  my $input = shift;
  my @result = sort {
    ($a->{seq} cmp $b->{seq})
    || ($a->{start} <=> $b->{start})
    || ($a->{end}   <=> $b->{end})
    || ($a->{score} <=> $b->{score})
  } @$input;
  return \@result;
}

sub by_seqandorientandpos_local {
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

sub by_modelandseqandorientandpos_local {
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

sub by_evalue_local {
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


1;
