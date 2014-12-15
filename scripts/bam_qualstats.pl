#!/usr/bin/perl
# Last changed Time-stamp: <2014-12-15 13:13:39 fabian>
#usage: perl bam_qualstats.pl -verbose -dir  ~/Work/ViennaNGS/Data/sinlge-end/ -odir /home/mescalin/fabian/Work/ViennaNGS/Progs/OuT


use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;
use Bio::Perl;
use Bio::DB::Sam;
use Statistics::R;

###########################
## # #   Variables   # # ##
###########################

my $bam_dir  = '';           # keeps path to directory with bams
my $odir     = '';           # keeps path to dir where results are stored
my $rlibpath = '/usr/bin/R'; # path to R installation (`which R`)
my %data     = ();           # stores all results from single bam files
my $VERBOSE  = 0;

my $match_control = 1;   # provids stats how many mapped bases match the reference genome
my $clip_control  = 1;   # provids stats how many bases are soft or hard clipped
my $split_control = 1;   # provids stats how many mapped reads are splitted (and how often); condideres CIGAR string, not good for segemehl)
my $qual_control  = 1;   # provids stats on quality of the match
my $edit_control  = 1;   # provids stats on the edit distance between read and mapped reference position
my $flag_control  = 1;   # analyses the sam bit flag for quality/strandedness/pair_vs_single end reads (must be =1)
my %raw_flag_data = ();
my $score_control = 0;   # provids stats on per-base quality scores
my $uniq_control  = 1;   # gives number and stats of multiplicity of readaligments (must be =1)






##################################
## # # Command-line Options # # ##
##################################

pod2usage(-verbose => 0)
        unless GetOptions(
	  "dir|d=s"     => \$bam_dir,
	  "odir|o=s"    => \$odir,
	  "rlib|r=s"    => \$rlibpath,
	  "match!"      => \$match_control,
	  "clip!"       => \$clip_control,
	  "qual!"       => \$qual_control,
	  "edit!"       => \$edit_control,
	  "help|h"      => sub{pod2usage(-verbose => 1)},
	  "man|m"       => sub{pod2usage(-verbose => 2)},
	  "verbose"     => sub{ $VERBOSE++ }
    );


## check/clean parameter and options
$bam_dir     =~ s/ /_/g;
$odir        =~ s/ /_/g;

unless ( -d "$bam_dir" ) {
  print STDERR "Directory $bam_dir not found\nEXIT\n";
  exit;
}

unless ( -e "$rlibpath" ) {
  print STDERR "Path to R ($rlibpath) not found\nEXIT\n";
  exit;
}

if ( -d "$odir" ) {
  my $expand_odir=`readlink -f $odir`;
  chomp($expand_odir);
  print STDERR "$expand_odir is used to store results. Existing files are overwritten.\n";
} else {
  print STDERR "Directory $odir not found\nEXIT\n";
  exit;
}

print STDERR "\nFollowing feature/attributes will be analysed:\n";
print STDERR "\t<sam flag>\t(paired- and single-end read counts; strand distribution)\n" if ($flag_control);
print STDERR "\t<uniqueness>\t(NH:i attribute; otherwise uniqueness is set to 1)\n" if ($uniq_control); 
print STDERR "\t<match>\t\t(distribution of percentage of read which matches refernce)\n" if ($match_control);
print STDERR "\t<clip>\t\t(distribution of percentage of read which are clipped during mapping)\n" if ($clip_control);
print STDERR "\t<edit>\t\t(distribution of edit distance between read and reference)\n" if ($edit_control);
print STDERR "\n";



##########################
## # #     MAIN     # # ##
##########################

##########################
## globs all bam files in $bam_dir
opendir(my $dh, $bam_dir) || die "can't opendir $bam_dir: $!";
my @bams = grep { /\.bam$/ && -f "$bam_dir/$_" } readdir($dh);

##########################
## make stat of single bam files
foreach my $bam (@bams){
  %{$data{$bam}}=&qual_singleBam("${bam_dir}/${bam}");
}

##########################
## print master results table
if ($uniq_control && $flag_control){
  print join("\t", 'sample', 'total alignments',  'mapped paires', 'unmapped paires', 'mapped singles', 'fwd strand',  'rev strand', 'mapped reads', 'uniq mapped reads')."\n";
  foreach my $sample (sort keys %data){
    print join("\t",
	       $sample,
	       $data{$sample}->{'aln_count'}->{'total'},
	       $data{$sample}->{'aln_count'}->{'mapped_pair'},
	       $data{$sample}->{'aln_count'}->{'unmapped_pair'},
	       $data{$sample}->{'aln_count'}->{'mapped_single'},
	       $data{$sample}->{'strand'}->{'forward'},
	       $data{$sample}->{'strand'}->{'reverse'},
	       $data{$sample}->{'uniq'}->{'mapped_reads'},
	       $data{$sample}->{'uniq'}->{'uniq_mapped_reads'})."\n";
  }
}

##########################
## plot boxplot for match
if ($match_control){
  my @Rstat_values_match = ();
  my @Rstat_names_match  = ();
  my @Rstat_ns_match     = ();
  my $Rstat_mcol_match   = 0;
  
  foreach my $sample (sort keys %data) {
    push @Rstat_values_match, $data{$sample}->{'match'}->{'min'};
    push @Rstat_values_match, $data{$sample}->{'match'}->{'1q'};
    push @Rstat_values_match, $data{$sample}->{'match'}->{'med'};
    push @Rstat_values_match, $data{$sample}->{'match'}->{'3q'};
    push @Rstat_values_match, $data{$sample}->{'match'}->{'max'};
    push @Rstat_names_match,  "'$sample'";
    push @Rstat_ns_match,     $data{$sample}->{'aln_count'}->{'total'};
    $Rstat_mcol_match++;
  }

  my $data_string_match = "list(stats=matrix(c(".join(",",@Rstat_values_match)."),5,$Rstat_mcol_match), n=c(".join(",",@Rstat_ns_match)."), names=c(".join(",",@Rstat_names_match)."))";
  &plot_boxplot("${odir}/match_stats", "match [%]", $data_string_match);
}

###########################
## plot barplot for counts
if ($flag_control){
  my @Rstat_data_count = ();
  
  # collect data for read.table string
  push @Rstat_data_count, 'samples', sort keys %data;                       # first line with sample names
  $Rstat_data_count[-1]="$Rstat_data_count[-1]\n";                          # ad new line to last column entry
  push @Rstat_data_count, 'mapped_single';               
  foreach my $sample (sort keys %data) {                                    # second line with all mapped_single count values
    push @Rstat_data_count, $data{$sample}->{'aln_count'}->{'mapped_single'};
  }
  $Rstat_data_count[-1]="$Rstat_data_count[-1]\n";                          # ad new line to last column entry
  push @Rstat_data_count, 'mapped_pair';
  foreach my $sample (sort keys %data) {                                    # third line with all mapped_pair count values
    push @Rstat_data_count, $data{$sample}->{'aln_count'}->{'mapped_pair'};
  }
  $Rstat_data_count[-1]="$Rstat_data_count[-1]\n";                          # ad new line to last column entry
  push @Rstat_data_count, 'unmapped_pair';
  foreach my $sample (sort keys %data) {                                    # fourth line with all unmapped_pair count values
    push @Rstat_data_count, $data{$sample}->{'aln_count'}->{'unmapped_pair'};
  }
  $Rstat_data_count[-1]="$Rstat_data_count[-1]\n";                          # ad new line to last column entry
  my $data_string_count = join(" ", @Rstat_data_count);                     # make read.table string
  &plot_barplot("${odir}/count_stats", "mapped reads", $data_string_count); # produce plot with read.table string input
}

##########################
## plot boxplot hard clip
if($clip_control){
  my @Rstat_values_clipH = ();
  my @Rstat_names_clipH  = ();
  my @Rstat_ns_clipH     = ();
  my $Rstat_mcol_clipH   = 0;
  
  foreach my $sample (sort keys %data) {
    push @Rstat_values_clipH, $data{$sample}->{'clip'}->{'hardclip'}->{'min'};
    push @Rstat_values_clipH, $data{$sample}->{'clip'}->{'hardclip'}->{'1q'};
    push @Rstat_values_clipH, $data{$sample}->{'clip'}->{'hardclip'}->{'med'};
    push @Rstat_values_clipH, $data{$sample}->{'clip'}->{'hardclip'}->{'3q'};
    push @Rstat_values_clipH, $data{$sample}->{'clip'}->{'hardclip'}->{'max'};
    push @Rstat_names_clipH,  "'$sample'";
    push @Rstat_ns_clipH,     $data{$sample}->{'aln_count'}->{'total'};
    $Rstat_mcol_clipH++;
  }

  my $data_string_clipH = "list(stats=matrix(c(".join(",",@Rstat_values_clipH)."),5,$Rstat_mcol_clipH), n=c(".join(",",@Rstat_ns_clipH)."), names=c(".join(",",@Rstat_names_clipH)."))";
  &plot_boxplot("${odir}/hardclipped_stats", "hard clipped [nt]", $data_string_clipH);
}

##########################
## plot boxplot soft clip
if($clip_control){
  my @Rstat_values_clipS = ();
  my @Rstat_names_clipS  = ();
  my @Rstat_ns_clipS     = ();
  my $Rstat_mcol_clipS   = 0;
  
  foreach my $sample (sort keys %data) {
    push @Rstat_values_clipS, $data{$sample}->{'clip'}->{'softclip'}->{'min'};
    push @Rstat_values_clipS, $data{$sample}->{'clip'}->{'softclip'}->{'1q'};
    push @Rstat_values_clipS, $data{$sample}->{'clip'}->{'softclip'}->{'med'};
    push @Rstat_values_clipS, $data{$sample}->{'clip'}->{'softclip'}->{'3q'};
    push @Rstat_values_clipS, $data{$sample}->{'clip'}->{'softclip'}->{'max'};
    push @Rstat_names_clipS,  "'$sample'";
    push @Rstat_ns_clipS,     $data{$sample}->{'aln_count'}->{'total'};
    $Rstat_mcol_clipS++;
  }

  my $data_string_clipS = "list(stats=matrix(c(".join(",",@Rstat_values_clipS)."),5,$Rstat_mcol_clipS), n=c(".join(",",@Rstat_ns_clipS)."), names=c(".join(",",@Rstat_names_clipS)."))";
  &plot_boxplot("${odir}/softclipped_stats", "soft clipped [nt]", $data_string_clipS);
}

##########################
## plot boxplot edit distance
if($edit_control){
  my @Rstat_values_edit = ();
  my @Rstat_names_edit  = ();
  my @Rstat_ns_edit     = ();
  my $Rstat_mcol_edit   = 0;
  
  foreach my $sample (sort keys %data) {
    push @Rstat_values_edit, $data{$sample}->{'edit'}->{'min'};
    push @Rstat_values_edit, $data{$sample}->{'edit'}->{'1q'};
    push @Rstat_values_edit, $data{$sample}->{'edit'}->{'med'};
    push @Rstat_values_edit, $data{$sample}->{'edit'}->{'3q'};
    push @Rstat_values_edit, $data{$sample}->{'edit'}->{'max'};
    push @Rstat_names_edit,  "'$sample'";
    push @Rstat_ns_edit,     $data{$sample}->{'aln_count'}->{'total'};
    $Rstat_mcol_edit++;
  }

  my $data_string_edit = "list(stats=matrix(c(".join(",",@Rstat_values_edit)."),5,$Rstat_mcol_edit), n=c(".join(",",@Rstat_ns_edit)."), names=c(".join(",",@Rstat_names_edit)."))";
  &plot_boxplot("${odir}/editdistance_stats", "edit distance", $data_string_edit);
}

##########################
## plot boxplot quality score
if($qual_control){
  my @Rstat_values_qual = ();
  my @Rstat_names_qual  = ();
  my @Rstat_ns_qual     = ();
  my $Rstat_mcol_qual   = 0;
  
  foreach my $sample (sort keys %data) {
    push @Rstat_values_qual, $data{$sample}->{'qual'}->{'min'};
    push @Rstat_values_qual, $data{$sample}->{'qual'}->{'1q'};
    push @Rstat_values_qual, $data{$sample}->{'qual'}->{'med'};
    push @Rstat_values_qual, $data{$sample}->{'qual'}->{'3q'};
    push @Rstat_values_qual, $data{$sample}->{'qual'}->{'max'};
    push @Rstat_names_qual,  "'$sample'";
    push @Rstat_ns_qual,     $data{$sample}->{'aln_count'}->{'total'};
    $Rstat_mcol_qual++;
  }

  my $data_string_qual = "list(stats=matrix(c(".join(",",@Rstat_values_qual)."),5,$Rstat_mcol_qual), n=c(".join(",",@Rstat_ns_qual)."), names=c(".join(",",@Rstat_names_qual)."))";
  &plot_boxplot("${odir}/qualityscore_stats", "read quality score", $data_string_qual);
}

#print Dumper(%data) if($VERBOSE);



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

  my @match_data    = ();  
  my %clip_data     = ();  
  my %split_data    = ();  
  my @qual_data     = ();  
  my @edit_data     = ();
  my %flag_data     = ();
  my @score_data    = ();
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
    $data_out{'aln_count'}->{'mapped_single'}  = $flag_data{'paired'}->{'single-end'};
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
      my $split_counts=( defined($split_data{'splits'}) )?($split_data{'splits'}):(0);
      $data_out{'split'}->{'distribution_of_multisplit'}->{"$splits"}=$split_counts;
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

sub plot_boxplot { #plot boxplot from min, 1q, median, 3q, max
  my $filename    = shift;
  my $ylab        = shift;
  my $data_string = shift;

  my $R = Statistics::R->new();
  $R->startR;
  $R->set('rlib', $rlibpath);
  $R->set('log_dir', $odir);
  $R->run("postscript('${filename}.eps')") ;
  $R->run("summarydata<-$data_string") ;
  $R->run("bxp(summarydata, medcol = 'red', ylab='$ylab', xlab='samples')") ;
  $R->run("dev.off()") ;
  $R->stopR;
}

sub plot_barplot { #plot barplot read.table text string
  my $filename    = shift;
  my $ylab        = shift;
  my $data_string = shift;

  my $R = Statistics::R->new();
  $R->startR;
  $R->set('rlib', $rlibpath);
  $R->set('log_dir', $odir);
  $R->run("postscript('${filename}.eps')") ;
  $R->run("dat<-read.table(text = \"$data_string\", header = TRUE, row.names=1)") ;
  $R->run("dat_m<-as.matrix(dat)") ;
  $R->run("barplot(dat_m, xlim=c(0,ncol(dat_m)+3), col=1:nrow(dat_m), legend.text = TRUE, args.legend = list(x = ncol(dat_m) + 3, y=max(colSums(dat_m)), bty = 'n' ), ylab='$ylab', xlab='samples')") ;
  $R->run("dev.off()") ;
  $R->stopR;
}

###########################
## # #   Variables   # # ##
###########################
__END__

=head1 NAME

bam_qualstats.pl -- Generates quality statistics from specified BAM input files.

=head1 SYNOPSIS
bam_qualstats.pl --dir <PATH> --odir <PATH> [--match] [--clip] [--qual] [--edit] [--rlib <PATH>] [[--help]

=head1 OPTIONS

=over
    
=item B<--dir>

Path to directory with BAM files. All BAM files in this directory will be processed.
    
=item B<--odir>

Path to output directory. In this directory several output files will be created. Existing files with identiacal names will be over written.

=item B<--match>

Provids stats how many mapped bases match the reference genome.

=item B<--clip>

Provids stats how many bases are soft or hard clipped.

=item B<--qual>

Provids stats on quality of the read match against the reference sequence.

=item B<--edit>

Provids stats on the edit distance between read and mapped reference position.
    
=item B<--rlib -r>

Path to the R library.

=item B<--help -h>

Print short help.

=item B<--man>

Prints the manual page and exits.

=item B<--verbose>

Prints additional output.

=back

=head1 DESCRIPTION

This program screens a given directory for bam files and returns veriouse statistics for the found samples.

=head1 AUTHOR

Fabian Amman E<lt>fabian@tbi.univie.ac.atE<gt>

=cut
