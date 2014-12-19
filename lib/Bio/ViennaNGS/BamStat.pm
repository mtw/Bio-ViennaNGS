# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-20 00:31:16 mtw>

package Bio::ViennaNGS::BamStat;

use 5.12.0;
use version; our $VERSION = qv('0.12_07');
use Bio::DB::Sam 1.39;
use Moose;
use Carp;
use Data::Dumper;
use namespace::autoclean;

has 'bam' => (
	      is => 'rw',
	      isa => 'Str',
	      required => 1,
	      predicate => 'has_bam',
	     );

has 'data' => (
	       is => 'ro',
	       isa => 'HashRef',
	       predicate => 'has_data',
	      );

has 'control_match' => ( # provides stats how many mapped bases match the reference genome
			is => 'rw',
			isa => 'Bool',
			default => '1',
			clearer => 'clear_control_match',
			predicate => 'has_control_match',
		       );

has 'control_clip' => ( # provides stats how many bases are soft or hard clipped
		       is => 'rw',
		       isa => 'Bool',
		       default => '1',
		       clearer => 'clear_control_clip',
		       predicate => 'has_control_clip',
		       );

has 'control_split' => ( # provides stats how many/often  mapped reads are split
			is => 'rw',
			isa => 'Bool',
			default => '1',
			clearer => 'clear_control_split',
			predicate => 'has_control_split',
		       );

has 'control_qual' => ( # provides stats on quality of the match
		       is => 'rw',
		       isa => 'Bool',
		       default => '1',
		       clearer => 'clear_control_qual',
		       predicate => 'has_control_qual',
		      );

has 'control_edit' => ( # provides stats on the edit distance between read and mapped reference
		       is => 'rw',
		       isa => 'Bool',
		       default => '1',
		        clearer => 'clear_control_edit',
		       predicate => 'has_control_edit',
		      );

has 'control_flag' => ( # analyses the sam bit flag for qual/strands/pair_vs_single reads
		       is => 'rw',
		       isa => 'Bool',
		       default => '1',
		       clearer => 'clear_control_flag',
		       predicate => 'has_control_flag',
		      );

has 'control_score' => ( # provides stats on per-base quality scores
			is => 'rw',
			isa => 'Bool',
			default => '1',
			clearer => 'clear_control_score',
			predicate => 'has_control_score',
		       );

has 'control_uniq' => ( #  gives number and stats of multiplicity of readaligments
		       is => 'rw',
		       isa => 'Bool',
		       default => '1',
		       clearer => 'clear_control_uniq',
		       predicate => 'has_control_uniq',
		      );

has 'is_segemehl' => ( # toggles to consider segemehl specific bam feature
		      is => 'rw',
		      isa => 'Bool',
		      default => '0',
		      predicate => 'has_is_segemehl',
		     );

has 'data_out' => (
		     is => 'rw',
		     isa => 'HashRef',
		     default => sub { {} },
		     predicate => 'has_data_out',
		    );

has 'data_match' => (
		     is => 'rw',
		     isa => 'ArrayRef',
		     default => sub { [] },
		     predicate => 'has_data_match',
		    );

has 'data_qual' => (
		    is => 'rw',
		    isa => 'ArrayRef',
		    default => sub { [] },
		    predicate => 'has_data_qual',
		   );

has 'data_edit' => (
		    is => 'rw',
		    isa => 'ArrayRef',
		    default => sub { [] },
		    predicate => 'has_data_edit',
		   );

has 'data_score' => (
		    is => 'rw',
		    isa => 'ArrayRef',
		    default => sub { [] },
		    predicate => 'has_data_score',
		   );

has 'data_clip' => (
		    is => 'rw',
		    isa => 'ArrayRef',
		    default => sub { [] },
		    predicate => 'has_data_clip',
		   );

has 'data_split' => (
		     is => 'rw',
		     isa => 'HashRef',
		     default => sub { {} },
		     predicate => 'has_data_split',
		    );

has 'data_flag' => (
		    is => 'rw',
		    isa => 'HashRef',
		    default => sub { {} },
		    predicate => 'has_data_flag',
		   );

has 'data_uniq' => (
		    is => 'rw',
		    isa => 'HashRef',
		    default => sub { {} },
		    predicate => 'has_data_uniq',
		   );

sub stat_singleBam {
  my ($self) = @_;
  my $this_function = (caller(0))[3];
  confess "ERROR [$this_function] Attribute 'bam' not found $!"
    unless ($self->has_bam);

  ## read in bam file using Bio::DB::Sam
  #my $sam = Bio::DB::Sam->new(-bam   => $self->bam);
  my $bam = Bio::DB::Bam->open($self->bam);

  my $header       = $bam->header;
  my $target_count = $header->n_targets;
  my @chromosomes = @{$header->target_name};
  my $target_names = \@chromosomes;

  while (my $align = $bam->read1) {
    my %flags = ();
    my $qname      = $align->qname;
    my $seqid      = $target_names->[$align->tid];
    my $start      = $align->pos+1;
    my $end        = $align->calend;
    my $strand     = $align->strand;
    my $cigar      = $align->cigar_str;
    my @scores     = $align->qscore;
    my $match_qual = $align->qual;
    my $flag       = $align->flag;
    my $attributes = $align->aux;
    #    print "\$seqid: $seqid\t \$qname: $qname\t $start: $start\t \$end: $end\t \$cigar: $cigar\t\$strand: $strand\t \$match_qual: $match_qual\t \$attributes: $attributes\n";
    my @tags = $align->get_all_tags;
    #    print ">> tags: @tags\n";

    #############################
    ### flag
    if ($self->has_control_flag){

      my $multisplit_weight = 1;
      if ($self->has_is_segemehl && $align->has_tag('XL') ) {
	$multisplit_weight = (1/($align->aux_get('XL')));
      }

      if( $align->get_tag_values('PAIRED')) {
	${$self->data_flag}{'paired'}->{'paired-end'}+=$multisplit_weight;
      }else{
	${$self->data_flag}{'paired'}->{'single-end'}+=$multisplit_weight;
      }

      if( $align->get_tag_values('REVERSED')) {
	${$self->data_flag}{'strand'}->{'reverse'}+=$multisplit_weight;
      }else{
	${$self->data_flag}{'strand'}->{'forward'}+=$multisplit_weight;
      }

      if( $align->get_tag_values('PAIRED') && $align->get_tag_values('M_UNMAPPED')){
	${$self->data_flag}{'pairs'}->{'unmapped_pair'}+=$multisplit_weight;
      }
      elsif( $align->get_tag_values('PAIRED') && !$align->get_tag_values('M_UNMAPPED')) {
	${$self->data_flag}{'pairs'}->{'mapped_pair'}+=$multisplit_weight;
      }

      if( $align->get_tag_values('DUPLICATE') ){
	${$self->data_flag}{'unmapped'}->{'duplicated'}+=$multisplit_weight;
      }
      if( $align->get_tag_values('QC_FAILED') ){
	${$self->data_flag}{'unmapped'}->{'qual_failed'}+=$multisplit_weight;
      }
      if( $align->get_tag_values('UNMAPPED') ){
	${$self->data_flag}{'unmapped'}->{'unmapped'}+=$multisplit_weight;
      }
    }

    #############################
    ### uniq
    if ($self->has_control_uniq){

      my $multisplit_weight = 1;
      if ($self->has_is_segemehl && $align->has_tag('XL') ) {
	$multisplit_weight = ($align->aux_get('XL'));
      }

      my $multimap_weight = 1;
      if ( $align->has_tag('NH') ) {
	$multimap_weight = ($align->aux_get('NH'));
      }

      ${$self->data_uniq}{'uniq'}+=(1/$multisplit_weight) if($multimap_weight == 1);
      ${$self->data_uniq}{'mapped'}+=1/($multimap_weight*$multisplit_weight);
      ${$self->data_uniq}{'multiplicity'}->{$multimap_weight}+=(1/$multisplit_weight);

    }

    #############################
    ### Quality Score read
    if($self->has_control_qual){
      if($self->has_is_segemehl && $match_qual == 255){
	#carp "WARN [$this_function] no match_qual for segemehl";
	$self->clear_control_qual;
      }
      elsif($match_qual){
	push @{$self->data_qual}, $match_qual;
      }
      else{
	carp "WARN [$this_function] No matchqual for read $qname";
	carp "WARN [$this_function] Setting \$self->has_control_qual to zero";
	$self->clear_control_qual;
      }
    }

    #############################
    #### Edit distance
    if($self->has_control_edit){
      if( $align->has_tag('NM') ){
	push @{$self->data_edit}, $align->aux_get('NM');
      }
      else{
	carp "WARN [$this_function] No <NM> attribute for read $qname";
	carp "WARN [$this_function] Setting \$self->has_control_edit to zero";
	$self->clear_control_edit;
      }
    }

    #############################
    ### Match bases
    if($self->has_control_match){
      if($cigar){
	push @{$self->data_match},
	  sprintf("%.2f", 100*(&cigarmatch($cigar))/(&cigarlength($cigar)));
      }
      else{
	carp "WARN [$this_function] No CIGAR string for read $qname";
	carp "WARN [$this_function] Setting \$self->has_control_match to zero";
	$self->clear_control_match;
      }
    }

    #############################
    ### Clip reads
    if($self->control_clip){
      if($cigar){
	push @{$self->data_clip}, sprintf("%.2f", 100*(&cigarclip($cigar)));
      }
      else{
	carp "WARN [$this_function] No CIGAR string for read $qname";
	carp "WARN [$this_function] Setting \$self->has_control_clip to zero";
	$self->clear_control_clip;
        }
      }

    #############################
    ### Split reads
    if($self->control_split){
      if ($self->has_is_segemehl && $align->has_tag('XL') ) {
	my $split_counts = $align->aux_get('XL');
	${$self->data_split}{$split_counts}+=1/($split_counts);
	${$self->data_split}{'total'}+=1/($split_counts);
      }
      elsif ( $self->has_is_segemehl ){ }
      else{
	if($cigar){
	  my $split_counts = cigarsplit($cigar);
	  ${$self->data_split}{$split_counts}++;
	  ${$self->split_data}{'total'}++;
	}
	else{
	  carp "WARN [$this_function] No CIGAR string for read $qname";
	  carp "WARN [$this_function] Setting \$self->has_control_split to zero";
	  $self->clear_control_split;
	}
      }
    }

  } # end while


  #############################
  ### Extract reference genome/chromosome sizes from BAM header
  for (my $i=0;$i<$header->n_targets;$i++){
    my $seqid = $header->target_name->[$i];
    my $length = $header->target_len->[$i];
    ${$self->data_out}{'chrom'}->{$seqid}=$length;
 }

  #############################
  #### uniq
  if ($self->has_control_uniq){
    ${$self->data_out}{'uniq'}->{'uniq_mapped_reads'}=${$self->data_uniq}{'uniq'};
    ${$self->data_out}{'uniq'}->{'mapped_reads'}=sprintf("%d", ${$self->data_uniq}{'mapped'});
    ${$self->data_out}{'uniq'}->{'uniq_mapped_reads_percent'}=
      sprintf("%.2f", (100*${$self->data_uniq}{'uniq'}/${$self->data_uniq}{'mapped'}));
    foreach my $multiplicity (sort {$a cmp $b} keys %{${$self->data_uniq}{'multiplicity'}}) {
      ${$self->data_out}{'uniq'}->{'distribution_of_multimapper'}->{"<$multiplicity>"}=
	(${$self->data_uniq}{'multiplicity'}->{$multiplicity}/$multiplicity);
    }
  }

  #############################
  ##### quality stats
  if($self->has_control_qual){
    my %qual_stats=%{stats(@${$self->data_qual})};
    if ($qual_stats{'min'} == 255 && $qual_stats{'max'} == 255){
      carp "WARN [$this_function] No matchqual (all values set to 255)";
      carp "WARN [$this_function] Setting \$self->has_control_qual to zero";
      $self->clear_control_qual;
    }
    else{
      ${$self->data_out}{'qual'}->{'min'}  = $qual_stats{'min'};
      ${$self->data_out}{'qual'}->{'1q'}   = $qual_stats{'1q'};
      ${$self->data_out}{'qual'}->{'mean'} = $qual_stats{'mean'};
      ${$self->data_out}{'qual'}->{'med'}  = $qual_stats{'med'};
      ${$self->data_out}{'qual'}->{'3q'}   = $qual_stats{'3q'};
      ${$self->data_out}{'qual'}->{'max'}  = $qual_stats{'max'};
    }
  }

  #############################
  ### Clip data
  if($self->has_control_clip){
    my %clip_stats=%{stats(@{$self->data_clip})};

    ${$self->data_out}{'clip'}->{'min'}  = $clip_stats{'min'};
    ${$self->data_out}{'clip'}->{'1q'}   = $clip_stats{'1q'};
    ${$self->data_out}{'clip'}->{'mean'} = $clip_stats{'mean'};
    ${$self->data_out}{'clip'}->{'med'}  = $clip_stats{'med'};
    ${$self->data_out}{'clip'}->{'3q'}   = $clip_stats{'3q'};
    ${$self->data_out}{'clip'}->{'max'}  = $clip_stats{'max'};
  }

  #############################
  ##### Match data
  if($self->has_control_match){
    my %match_stats=%{stats(@{$self->data_match})};

    ${$self->data_out}{'match'}->{'min'}  = $match_stats{'min'};
    ${$self->data_out}{'match'}->{'1q'}   = $match_stats{'1q'};
    ${$self->data_out}{'match'}->{'mean'} = $match_stats{'mean'};
    ${$self->data_out}{'match'}->{'med'}  = $match_stats{'med'};
    ${$self->data_out}{'match'}->{'3q'}   = $match_stats{'3q'};
    ${$self->data_out}{'match'}->{'max'}  = $match_stats{'max'};
  }

  #############################
  ##### Edit distance
  if($self->has_control_edit){
    my %edit_stats=%{stats(@{$self->data_edit})};

    ${$self->data_out}{'edit'}->{'min'}  = $edit_stats{'min'};
    ${$self->data_out}{'edit'}->{'1q'}   = $edit_stats{'1q'};
    ${$self->data_out}{'edit'}->{'mean'} = $edit_stats{'mean'};
    ${$self->data_out}{'edit'}->{'med'}  = $edit_stats{'med'};
    ${$self->data_out}{'edit'}->{'3q'}   = $edit_stats{'3q'};
    ${$self->data_out}{'edit'}->{'max'}  = $edit_stats{'max'};
  }

  #############################
  ###### flags
  if ($self->has_control_flag){

    ${$self->data_flag}{'pairs'}->{'mapped_pair'}   = 0
      unless ( defined(${$self->data_flag}{'pairs'}->{'mapped_pair'})   );
    ${$self->data_flag}{'pairs'}->{'unmapped_pair'} = 0
      unless ( defined(${$self->data_flag}{'pairs'}->{'unmapped_pair'}) );

    ${$self->data_flag}{'paired'}->{'paired-end'}   = 0
      unless ( defined(${$self->data_flag}{'paired'}->{'paired-end'})   );
    ${$self->data_flag}{'paired'}->{'single-end'}   = 0
      unless ( defined(${$self->data_flag}{'paired'}->{'single-end'})   );

    my $total_mapped  = ${$self->data_flag}{'paired'}->{'paired-end'} +
      ${$self->data_flag}{'paired'}->{'single-end'};

    ${$self->data_out}{'aln_count'}->{'mapped_pair'} =
      ${$self->data_flag}{'pairs'}->{'mapped_pair'};
    ${$self->data_out}{'aln_count'}->{'unmapped_pair'} =
      ${$self->data_flag}{'pairs'}->{'unmapped_pair'};
    ${$self->data_out}{'aln_count'}->{'mapped_single'} =
      ${$self->data_flag}{'paired'}->{'single-end'};
    ${$self->data_out}{'aln_count'}->{'total'} = $total_mapped;

    ${$self->data_flag}{'strand'}->{'forward'}   = 0
      unless ( defined(${$self->data_flag}{'strand'}->{'forward'})   );
    ${$self->data_flag}{'strand'}->{'reverse'}   = 0 
      unless ( defined(${$self->data_flag}{'strand'}->{'reverse'})   );

    my $total_strands = ${$self->data_flag}{'strand'}->{'forward'} + 
      ${$self->data_flag}{'strand'}->{'reverse'};

    ${$self->data_out}{'strand'}->{'forward'} = ${$self->data_flag}{'strand'}->{'forward'};
    ${$self->data_out}{'strand'}->{'reverse'} = ${$self->data_flag}{'strand'}->{'reverse'};
  }

  #############################
  ###### split
  if ($self->has_control_split){
    foreach my $splits (sort {$a cmp $b} keys %{$self->data_split}) {
      my $split_counts=( defined(${$self->data_split}{$splits}) )?(${$self->data_split}{$splits}):(0);
      ${$self->data_out}{'split'}->{'distribution_of_multisplit'}->{"$splits"}=$split_counts;
    }
  }

}

__PACKAGE__->meta->make_immutable;


no Moose;

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

sub cigarclip { #usage: cigarclip($cigarstring)
  my $cigar_string = shift;
  my $cigar_length = 0;
  my $cigar_clip   = 0;
  while($cigar_string=~m/(\d+)[SH]/g) {$cigar_clip+=$1}
  while($cigar_string=~m/(\d+)\D/g){$cigar_length+=$1}
  return ($cigar_length==0)?(0):($cigar_clip/$cigar_length);
}



1;
__END__


=head1 NAME

Bio::ViennaNGS::BamStatSingle - Moose interface to BAM mapping statistics

=head1 SYNOPSIS

  use Bio::ViennaNGS::BamStatSingle;

  my $bss1 = Bio::ViennaNGS::BamStatSingle->new(bam => "path/to/file.bam");

=head1 DESCRIPTION

This module provides a L<Moose> interface to the mapping statistics of
a single BAM file. It builds on L<Bio::DB::Sam> and 

=head1 SEE ALSO

=over 

=item L<Bio::ViennaNGS>

=back

=head1 AUTHORS

=over 

=item Fabian Amman E<lt>fabian@tbi.univie.ac.atE<gt>

=item Michael T. Wolfinger  E<lt>michael@wolfinger.euE<gt>

=back

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014 by Michael T. Wolfinger

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.16.3 or,
at your option, any later version of Perl 5 you may have available.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut
