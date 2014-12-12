#!/usr/bin/perl
# Last changed Time-stamp: <2014-12-12 15:02:14 fabian>



# Extracts quality information from bam/sam files
# Usage ./quality.pl <bam.file>


use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Bio::Perl;
use Bio::DB::Sam;






##########################
## # #     MAIN     # # ##
##########################

my $bam_dir = shift;
opendir(my $dh, $bam_dir) || die "can't opendir $bam_dir: $!";
my @bams = grep { /\.bam$/ && -f "$bam_dir/$_" } readdir($dh);

my %data = ();

print join("\t", @bams)."\n";
foreach my $bam (@bams){
  print "\n################\n## $bam\n\n";
  %{$data{$bam}}=&qual_singleBam("${bam_dir}/${bam}");
}

print Dumper(%data);

###########################
## # #  subroutines  # # ##
###########################



sub qual_singleBam{

  ## data structure to return from subroutine
  my %data_out  = ();
  
  ## get input file
  my $bam_fh   = shift;

  ## sets variables and contols if 
  ## feature should be considered

  my $match_control = 1;   # provids stats how many mapped bases match the reference genome
  my @match_data    = ();
  
  my $clip_control  = 1;   # provids stats how many bases are soft or hard clipped
  my %clip_data     = ();
  
  my $split_control = 1;   # provids stats how many mapped reads are splitted (and how often); condideres CIGAR string, not good for segemehl)
  my %split_data    = ();
  
  my $qual_control  = 1;   # provids stats on quality of the match
  my @qual_data     = ();
  
  my $edit_control  = 1;   # provids stats on the edit distance between read and mapped reference position
  my @edit_data     = ();
  
  my $flag_control  = 1;   # analyses the sam bit flag for quality/strandedness/pair_vs_single end reads
  my %raw_flag_data = ();
  my %flag_data     = ();
  
  my $score_control = 1;   # provids stats on per-base quality scores
  my @score_data    = ();
  
  my $uniq_control  = 1;   # gives number and stats of multiplicity of readaligments
  my %uniq_data     = ('uniq', 0, 'mapped', 0);
  
  
  ## read in bam file using Bio::DB::Sam
  my $sam = Bio::DB::Sam->new(-bam   => "$bam_fh");
  
  
  ## get chromosomes name
  my @chromosomes    = $sam->seq_ids;
  
  
  ## get single alignments and loop over each single alignment
  my @alignments = $sam->features();
  for my $a (@alignments) {
    my $qname      = $a->qname;                     # read id
    my $seqid      = $a->seq_id;                    # reference seq id
    my $start      = $a->start;                     # start position of mapping
    my $end        = $a->end;                       # end position of mapping
    my $strand     = $a->strand;                    # strand of mapping
    my $cigar      = $a->cigar_str;                 # CIGAR string
    my @scores     = $a->qscore;                    # per-base quality scores
    my $match_qual = $a->qual;                      # quality of the match
    my $flag       = $a->flag;                      # bit wise flag
    my $type       = $a->type;                      # get type (match, not-matched??)
    my (@features) = $a->get_SeqFeatures;           # no clue what kind of feature 
    my %attributes = $a->attributes;                # reports optinal attributest from the last sam column
    
    
    if (keys %attributes){
      my %flags      = map {$_ => 1} split('\|', $attributes{'FLAGS'});
      my $paired     = ( defined($flags{'PAIRED'}) )?(1):(0);
      
      
      #############################
      ### flag
      if ($flag_control){
        if( defined($flag) ){
          $raw_flag_data{$flag}++;
	  $raw_flag_data{'total'}++;
	  
          if( defined($flags{'PAIRED'}) ){
	    $flag_data{'paired'}->{'paired-end'}++;
	  }else{
	    $flag_data{'paired'}->{'single-end'}++;
	  }
	  
	  if( defined($flags{'REVERSED'}) ){
            $flag_data{'strand'}->{'reverse'}++;
          }else{
            $flag_data{'strand'}->{'forward'}++;
          }
	  
	  if( defined($flags{'PAIRED'}) && defined($flags{'M_UNMAPPED'}) ){
	    $flag_data{'pairs'}->{'mapped_pair'}++;
	  }
	  elsif( defined($flags{'PAIRED'}) && !defined($flags{'M_UNMAPPED'}) ){
	    $flag_data{'pairs'}->{'unmapped_pair'}++;
          }
	  
	  if( defined($flags{'DUPLICATE'}) ){
	    $flag_data{'unmapped'}->{'duplicated'}++;
	  }
	  if( defined($flags{'QC_FAILED'}) ){
	    $flag_data{'unmapped'}->{'qual_failed'}++;
	  }
	  if( defined($flags{'UNMAPPED'}) ){
	    $flag_data{'unmapped'}->{'unmapped'}++;
	  }
	  
        }
      }
      
      #############################
      ### uniq
      if ($uniq_control){
        if( defined($attributes{'NH'}) ){
          $uniq_data{'uniq'}++ if($attributes{'NH'} == 1);
          $uniq_data{'mapped'}+=1/$attributes{'NH'};
	  $uniq_data{'multiplicity'}->{$attributes{'NH'}}++;
	  #$uniq_data{'multiplicity'}->{'total'}++;
        }
	else{
	  $uniq_data{'uniq'}++;
	  $uniq_data{'mapped'}++;
	  $uniq_data{'multiplicity'}->{'1'}++;
	}
      }
      
      #############################
      ### Quality Score read
      if($qual_control){
        if($match_qual){
	  push @qual_data, $match_qual;
	}
	else{
	  print STDERR "warnings: no matchqual for read $qname\nwarnings: Setting \$qual_control to zero\n";
	  $qual_control = 0;
	}
      }

      #############################
      #### Edit distance
      if($edit_control){
        if( defined($attributes{'NM'}) ){  
	  push @edit_data, $attributes{'NM'};
	}
        else{
          print STDERR "warnings: no <NM> attribute for read $qname\nwarnings: Setting \$edit_control to zero\n";
          $edit_control = 0;
        }
      }


      #############################
      ### Match bases
      if($match_control){
        if($cigar){
          push @match_data, sprintf("%.2g", 100*(&cigarmatch($cigar))/(&cigarlength($cigar)));
        }
        else{
          print STDERR "warnings: no CIGAR string for read $qname\nwarnings: Setting \$match_control to zero\n";
          $match_control = 0;
        }
      }

      #############################
      ### Clip reads
      if($clip_control){
        if($cigar){
          push @{$clip_data{'H'}}, sprintf("%.2g", 100*(&cigarHclip($cigar)));
          push @{$clip_data{'S'}}, sprintf("%.2g", 100*(&cigarSclip($cigar)));
        }
        else{
          print STDERR "warnings: no CIGAR string for read $qname\nwarnings: Setting \$clip_control to zero\n";
          $clip_control = 0;
        }
      }
 
  
      #############################
      ### Split reads
      if($split_control){
        if($cigar){
          $split_data{&cigarsplit($cigar)}++;
          $split_data{'total'}++;
        }
        else{
          print STDERR "warnings: no CIGAR string for read $qname\nwarnings: Setting \$split_control to zero\n";
          $split_control = 0;
        }
      }
      
    }
  }
  
  foreach my $chrom (sort @chromosomes){
    my $length = $sam->length($chrom);
    $data_out{'chrom'}->{$chrom}=$length;
  }
  
  #############################
  #### uniq
  if ($uniq_control){
    $data_out{'uniq'}->{'uniq_mapped_reads'}=$uniq_data{'uniq'};
    $data_out{'uniq'}->{'mapped_reads'}=sprintf("%d", $uniq_data{'mapped'});
    $data_out{'uniq'}->{'uniq_mapped_reads_percent'}=sprintf("%.2g", (100*$uniq_data{'uniq'}/$uniq_data{'mapped'}));
    foreach my $multiplicity (sort {$a cmp $b} keys %{$uniq_data{'multiplicity'}}) {
      $data_out{'uniq'}->{'distribution_of_multimapper'}->{"<$multiplicity>"}=($uniq_data{'multiplicity'}->{$multiplicity}/$multiplicity);
    }
    
  }
  
  #############################
  ##### quality stats
  if($qual_control){
    my %qual_stats=%{&stats(@qual_data)};
    if ($qual_stats{'min'} == 255 && $qual_stats{'max'} == 255){
      print STDERR "warnings: no matchqual (all values set to 255)\nwarnings: Setting \$qual_control to zero\n";
      $qual_control = 0;
    }
    else{
      $data_out{'qual'}->{'min'}  = $qual_stats{'min'};
      $data_out{'qual'}->{'1q'}   = $qual_stats{'1q'};
      $data_out{'qual'}->{'mean'} = $qual_stats{'mean'};
      $data_out{'qual'}->{'med'}  = $qual_stats{'med'};
      $data_out{'qual'}->{'3q'}   = $qual_stats{'3q'};
      $data_out{'qual'}->{'max'}  = $qual_stats{'max'};
    }
  }
  
  #############################
  ### Clip data
  if($clip_control){
    my %clipH_stats=%{&stats(@{$clip_data{'H'}})};
    my %clipS_stats=%{&stats(@{$clip_data{'S'}})};

    $data_out{'clip'}->{'hardclip'}->{'min'}  = $clipH_stats{'min'};
    $data_out{'clip'}->{'hardclip'}->{'1q'}   = $clipH_stats{'1q'};
    $data_out{'clip'}->{'hardclip'}->{'mean'} = $clipH_stats{'mean'};
    $data_out{'clip'}->{'hardclip'}->{'med'}  = $clipH_stats{'med'};
    $data_out{'clip'}->{'hardclip'}->{'3q'}   = $clipH_stats{'3q'};
    $data_out{'clip'}->{'hardclip'}->{'max'}  = $clipH_stats{'max'};

    $data_out{'clip'}->{'softclip'}->{'min'}  = $clipS_stats{'min'};
    $data_out{'clip'}->{'softclip'}->{'1q'}   = $clipS_stats{'1q'};
    $data_out{'clip'}->{'softclip'}->{'mean'} = $clipS_stats{'mean'};
    $data_out{'clip'}->{'softclip'}->{'med'}  = $clipS_stats{'med'};
    $data_out{'clip'}->{'softclip'}->{'3q'}   = $clipS_stats{'3q'};
    $data_out{'clip'}->{'softclip'}->{'max'}  = $clipS_stats{'max'};

  }
  
  
  #############################
  ##### Match data
  if($match_control){
    my %match_stats=%{&stats(@match_data)};

    $data_out{'match'}->{'min'}  = $match_stats{'min'};
    $data_out{'match'}->{'1q'}   = $match_stats{'1q'};
    $data_out{'match'}->{'mean'} = $match_stats{'mean'};
    $data_out{'match'}->{'med'}  = $match_stats{'med'};
    $data_out{'match'}->{'3q'}   = $match_stats{'3q'};
    $data_out{'match'}->{'max'}  = $match_stats{'max'};

  }
  
  
  #############################
  ##### Edit distance
  if($edit_control){
    my %edit_stats=%{&stats(@edit_data)};

    $data_out{'edit'}->{'min'}  = $edit_stats{'min'};
    $data_out{'edit'}->{'1q'}   = $edit_stats{'1q'};
    $data_out{'edit'}->{'mean'} = $edit_stats{'mean'};
    $data_out{'edit'}->{'med'}  = $edit_stats{'med'};
    $data_out{'edit'}->{'3q'}   = $edit_stats{'3q'};
    $data_out{'edit'}->{'max'}  = $edit_stats{'max'};
 
    
  }
  
  
  #############################
  ###### flags
  if ($flag_control){
        
    $flag_data{'pairs'}->{'mapped_pair'}   = 0 unless ( defined($flag_data{'pairs'}->{'mapped_pair'})   );
    $flag_data{'pairs'}->{'unmapped_pair'} = 0 unless ( defined($flag_data{'pairs'}->{'unmapped_pair'}) );
    
    $flag_data{'paired'}->{'paired-end'}   = 0 unless ( defined($flag_data{'paired'}->{'paired-end'})   );
    $flag_data{'paired'}->{'single-end'}   = 0 unless ( defined($flag_data{'paired'}->{'single-end'})   );
    my $total_mapped                       = $flag_data{'paired'}->{'paired-end'} + $flag_data{'paired'}->{'single-end'};
    
    $data_out{'aln_count'}->{'mapped_pair'}    = $flag_data{'pairs'}->{'mapped_pair'};
    $data_out{'aln_count'}->{'unmapped_pair'}  = $flag_data{'pairs'}->{'unmapped_pair'};
    $data_out{'aln_count'}->{'mapped_single'}  = $flag_data{'pairs'}->{'single-end'};
    $data_out{'aln_count'}->{'total'}          = $total_mapped;

    $flag_data{'strand'}->{'forward'}   = 0 unless ( defined($flag_data{'strand'}->{'forward'})   );
    $flag_data{'strand'}->{'reverse'}   = 0 unless ( defined($flag_data{'strand'}->{'reverse'})   );
    my $total_strands                   = $flag_data{'strand'}->{'forward'} + $flag_data{'strand'}->{'reverse'};

    $data_out{'strand'}->{'forward'}    = $flag_data{'strand'}->{'forward'};
    $data_out{'strand'}->{'reverse'}    = $flag_data{'strand'}->{'reverse'};
  }
  
  
  #############################
  ###### split
  if ($split_control){
    foreach my $splits (sort {$a cmp $b} keys %split_data) {
      $data_out{'split'}->{'distribution_of_multisplit'}->{"<$splits>"}=$split_data{'splits'};
    }
  }

  return( %data_out )
}

sub stats{
  # usage: %h = %{stats(@a)};
  my @vals = sort {$a <=> $b} @_;
  my %stats = ();
  my $median = '';
  
  if(@vals%2){$stats{'med'} = $vals[int(@vals/2)];}                      #odd median
  else{$stats{'med'} = ($vals[int(@vals/2)-1] + $vals[int(@vals/2)])/2;} #even median
  $stats{'med'}             = sprintf("%.2f", $stats{'med'});
  $stats{'mean'}            = sprintf("%.2f", &mean(\@vals));            ## mean
  $stats{'min'}             = sprintf("%.2f", &min(\@vals));             ## min
  $stats{'max'}             = sprintf("%.2f", &max(\@vals));             ## max
  $stats{'1q'}              = sprintf("%.2f", $vals[int(@vals/4)]);      ## 1.quartile
  $stats{'3q'}              = sprintf("%.2f", $vals[int((@vals*3)/4)]);  ## 3.quartile

  return(\%stats);
}

sub mean { # usage: $h = %{mean(\@a)};
  my ($arrayref) = @_;
  my $sum;
  foreach (@$arrayref) {$sum += $_}
  return $sum / @$arrayref;
}

sub max { # usage: $h = %{max(\@a)};
  my ($arrayref) = @_;
  my $max = $arrayref->[0];
  foreach (@$arrayref) {$max = $_ if $_ > $max}
  return $max;
}

sub min { # usage: $h = %{min(\@a)};
  my ($arrayref) = @_;
  my $min = $arrayref->[0];
  foreach (@$arrayref) {$min = $_ if $_ < $min}
  return $min;
}

sub cigarlength { #usage: &cigarlength($cigarstring)
  my $cigar_string = shift;
  my $cigar_length = 0;
  while($cigar_string=~m/(\d+)[MIX=]/g){$cigar_length+=$1}
  return($cigar_length);
}

sub cigarmatch { #usage: cigarmatch($cigarstring)
  my $cigar_string = shift;
  my $cigar_match  = 0;
  while($cigar_string=~m/(\d+)[M=]/g){$cigar_match+=$1}
  return($cigar_match);
}

sub cigarsplit { #usage: cigarsplit($cigarstring)
  my $cigar_string = shift;
  my $cigar_split  = 0;
  while($cigar_string=~m/N/g){$cigar_split+=1}
  return($cigar_split);
}

sub cigarSclip { #usage: cigarSclip($cigarstring)
  my $cigar_string = shift;
  my $cigar_length = 0;
  my $cigar_clip   = 0;
  while($cigar_string=~m/(\d+)S/g) {$cigar_clip+=$1}
  while($cigar_string=~m/(\d+)\D/g){$cigar_length+=$1}
  return ($cigar_length==0)?(0):($cigar_clip/$cigar_length);
}

sub cigarHclip { #usage: cigarHclip($cigarstring)
  my $cigar_string = shift;
  my $cigar_clip   = 0;
  my $cigar_length = 0;
  while($cigar_string=~m/(\d+)H/g) {$cigar_clip+=$1}
  while($cigar_string=~m/(\d+)\D/g){$cigar_length+=$1}
  return ($cigar_length==0)?(0):($cigar_clip/$cigar_length);
}

	
