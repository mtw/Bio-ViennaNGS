# -*-CPerl-*-
# Last changed Time-stamp: <2015-06-29 15:40:26 mtw>

package Bio::ViennaNGS::BamStatSummary;

use version; our $VERSION = qv('0.15');
use Moose;
use Carp;
use POSIX qw(floor);
use Statistics::R;
use Data::Dumper;
use Path::Class;
use namespace::autoclean;
use Bio::ViennaNGS::BamStat;
use Tie::Hash::Indexed;
use File::Basename;

has 'data' => (
	       is => 'ro',
	       isa => 'ArrayRef [Bio::ViennaNGS::BamStat]',
	       default => sub { [] },
	      );

has 'countStat' => (
		    is => 'rw',
		    isa => 'HashRef',
		    predicate => 'has_countStat',
		    default => sub { {} },
#		    auto_deref => '1',
		   );

has 'outpath' => (
		  is => 'rw',
		  isa => 'Str',
		  required => '1',
		 );

has 'rlib' => (
	       is => 'rw',
	       isa => 'Str',
	       required => '0',
	      );

has 'files' => (
		is => 'rw',
		isa => 'ArrayRef',
		required => 1,
		predicate => 'has_files',
	       );

has 'control_match' => ( # provides stats how many mapped bases match the reference genome
			is => 'rw',
			isa => 'Bool',
			default => '0',
			predicate => 'has_control_match',
		       );

has 'control_clip' => ( # provides stats how many bases are soft or hard clipped
		       is => 'rw',
		       isa => 'Bool',
		       default => '0',
		       predicate => 'has_control_clip',
		       );

has 'control_split' => ( # provides stats how many/often  mapped reads are split
			is => 'rw',
			isa => 'Bool',
			default => '0',
			predicate => 'has_control_split',
		       );

has 'control_qual' => ( # provides stats on quality of the match
		       is => 'rw',
		       isa => 'Bool',
		       default => '0',
		       predicate => 'has_control_qual',
		      );

has 'control_edit' => ( # provides stats on the edit distance between read and mapped reference
		       is => 'rw',
		       isa => 'Bool',
		       default => '0',
		       predicate => 'has_control_edit',
		      );

has 'control_flag' => ( # analyses the sam bit flag for qual/strands/pair_vs_single reads
		       is => 'rw',
		       isa => 'Bool',
		       default => '0',
		       predicate => 'has_control_flag',
		      );

has 'control_score' => ( # provides stats on per-base quality scores
			is => 'rw',
			isa => 'Bool',
			default => '0',
			predicate => 'has_control_score',
		       );

has 'control_uniq' => ( #  gives number and stats of multiplicity of readaligments
		       is => 'rw',
		       isa => 'Bool',
		       default => '1',
		       predicate => 'has_control_uniq',
		      );

has 'is_segemehl' => ( # toggles to consider segemehl specific bam feature
		      is => 'rw',
		      isa => 'Bool',
		      default => '0',
		      predicate => 'has_is_segemehl',
		     );

sub populate_data {
  my ($self) = @_;

  foreach my $bamfile (@{$self->files}){
    #carp ">> processing $bamfile\n";
    
    my $bo = Bio::ViennaNGS::BamStat->new(bam            => $bamfile,
					  control_edit   => $self->control_edit,
					  control_flag   => $self->control_flag,
					  control_score  => $self->control_score,
					  control_match  => $self->control_match,
					  control_clip   => $self->control_clip,
					  control_split  => $self->control_split,
					  control_qual   => $self->control_qual,
					  control_uniq   => $self->control_uniq,
					  is_segemehl    => $self->is_segemehl,
					 );

    $bo->stat_singleBam();
    push (@{$self->data}, $bo);
  }
}


sub populate_countStat {
  my ($self) = @_;

  tie my %hdr, 'Tie::Hash::Indexed';
  %hdr = (
	  "sample"           => "Sample",
	  "total_alignments" => "# Total alignments",
	  "mapped_reads"     => "# Mapped reads",
	  "umapped_reads"    => "# Unique mapped reads",
	  "mmapped_reads"    => "# Multi mapped reads",
	  "aligned_pairs"    => "# Aligned in pairs",
	  "aligned_mm"       => "# Aligned mate missing",
	  "aligned_se"       => "# Aligned single end",
	  "aligned_fwd"      => "# Aligned forward strand",
	  "aligned_rev"      => "# Aligned reverse strand");
  ${$self->countStat}{'header'} = \%hdr;

  foreach my $sample (@{$self->data}){
    my ($basename,$dir,$ext) = fileparse($$sample{'bam'},qr/\.[^.]*/);
     ${$self->countStat}{$basename}{'total_alignments'} = floor(0.5 + $$sample{'data_out'}->{'aln_count'}->{'total'} );
     ${$self->countStat}{$basename}{'mapped_reads'}     = floor(0.5 + $$sample{'data_out'}->{'uniq'}->{'mapped_reads'} );
     ${$self->countStat}{$basename}{'umapped_reads'}    = floor(0.5 + $$sample{'data_out'}->{'uniq'}->{'uniq_mapped_reads'} );
     ${$self->countStat}{$basename}{'mmapped_reads'}    = floor(0.5 + $$sample{'data_out'}->{'uniq'}->{'mapped_reads'} - $$sample{'data_out'}->{'uniq'}->{'uniq_mapped_reads'} );
     ${$self->countStat}{$basename}{'aligned_pairs'}    = floor(0.5 + ($$sample{'data_out'}->{'aln_count'}->{'mapped_pair'})/2 );
     ${$self->countStat}{$basename}{'aligned_mm'}       = floor(0.5 + $$sample{'data_out'}->{'aln_count'}->{'unmapped_pair'} );
     ${$self->countStat}{$basename}{'aligned_se'}       = floor(0.5 + $$sample{'data_out'}->{'aln_count'}->{'mapped_single'} );
     ${$self->countStat}{$basename}{'aligned_fwd'}      = floor(0.5 + $$sample{'data_out'}->{'strand'}->{'forward'} );
     ${$self->countStat}{$basename}{'aligned_rev'}      = floor(0.5 + $$sample{'data_out'}->{'strand'}->{'reverse'} );
  }
}

sub dump_countStat {
  my ($self,$how) = @_;
  my $mn = "mapping_stats.csv";
  my $fn = file($self->outpath,$mn);
  
  open(OUT, "> $fn") or croak "cannot open OUT $!";
  print OUT join ("\t", values %{$self->countStat->{'header'} })."\n";

  foreach my $sample (keys %{$self->countStat} ){
    next if ($sample eq 'header');
    my @line = ();
    foreach my $key ( keys %{$self->countStat->{'header'}} ) {
      if ($key eq 'sample'){
	push (@line, $sample);
	next;
      }
      push @line, $self->countStat->{$sample}->{$key};
    }
    print OUT join ("\t", @line)."\n";
  }

  close (OUT);
}

sub make_BarPlot{
  my ($self) = @_;
  my @Rstat_data_count = ();

  ## collect data for read.table string
  push @Rstat_data_count, 'Samples', grep {!/header/} keys %{$self->countStat}; # first line with sample names
  $Rstat_data_count[-1]="$Rstat_data_count[-1]\n";                              # end 1st line

  push @Rstat_data_count, 'aligned_se';
  foreach my $sample ( grep {!/header/}  keys %{$self->countStat} ) {    # 2nd line with single end aln counts
    push @Rstat_data_count, $self->countStat->{$sample}{'aligned_se'};
  }
  $Rstat_data_count[-1]="$Rstat_data_count[-1]\n";                       # end 2nd line

  push @Rstat_data_count, 'aligned_pairs';
  foreach my $sample ( grep {!/header/}  keys %{$self->countStat} ) {    # 3rd line with aligned pairs count
    push @Rstat_data_count, $self->countStat->{$sample}{'aligned_pairs'};
  }
  $Rstat_data_count[-1]="$Rstat_data_count[-1]\n";                       # end 3rd line

  push @Rstat_data_count, 'aligned_mm';
  foreach my $sample ( grep {!/header/}  keys %{$self->countStat} ) {    # 4th line with incomplete aligned pairs
    push @Rstat_data_count, $self->countStat->{$sample}{'aligned_mm'};
  }
  $Rstat_data_count[-1]="$Rstat_data_count[-1]\n";                       # end 4th line

  ## produce bar plot
  my $mn = "mapping_stats.pdf";
  my $fn = file($self->outpath,$mn);
  my $datastring = join(" ", @Rstat_data_count);
  $self->plot_barplot($fn, "Mapped reads", $datastring); # produce plot with read.table string input
}

sub plot_barplot { #plot barplot read.table text string
  my ($self, $filename, $ylab, $data_string) = @_;
  my ($bn,$odir,$ext) = fileparse($filename, qr /\..*/);

  my $rlibpath        = $self->rlib;

  $filename .= '.pdf' unless ($ext eq '.pdf');

  my $R = Statistics::R->new();
  $R->startR;
  $R->set('rlib', $rlibpath);
  $R->set('log_dir', $odir);
  $R->run("pdf('${filename}')") ;
  $R->run("dat<-read.table(text = \"$data_string\", header = TRUE, row.names=1)") ;
  $R->run("dat_m<-as.matrix(dat)") ;
# $R->run("colors<-terrain.colors(nrow(dat_m), alpha = 1)") ;
  $R->run("colors<-c('lightblue','lightgreen','lightcoral', terrain.colors(nrow(dat_m)-3, alpha = 1))") ;
  $R->run("types<-row.names(dat_m)") ;
  $R->run("par(mar = c(15,3,5,5), oma = c(1, 1, 4, 1))") ;
# $R->run("barplot(dat_m, xlim=c(0,ncol(dat_m)+2), col=colors, legend.text = TRUE, args.legend = list(x = ncol(dat_m) + 2, y=max(colSums(dat_m)), bty = 'n' ), ylab='$ylab', xlab='Samples')") ;
# $R->run("barplot(dat_m, xlim=c(0,ncol(dat_m)), col=colors, legend.text = TRUE, args.legend = list(x = ncol(dat_m) + 5, y=-5, bty = 'o' ), ylab='$ylab', xlab='Samples', las=3)") ;
# $R->run("barplot(dat_m, xlim=c(0,ncol(dat_m)), col=colors, legend.text = TRUE, args.legend = list(\"topright\", horiz = TRUE, bty = 'o' ), ylab='$ylab', xlab='', las=3)") ;
  $R->run("barplot(dat_m, xlim=c(0,ncol(dat_m)), col=colors, ylab='$ylab', xlab='', las=3)") ;
  $R->run("par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0),mar = c(0, 0, 0, 0), new = TRUE)") ;
  $R->run("legend('top', types, horiz = TRUE, inset = c(0,0), bty = 'n', fill = colors, cex = 1.2 )") ;
  $R->run("dev.off()") ;
  $R->stopR;
}


sub make_BoxPlot{
  my ($self, $whattodo) = @_;
  my @Rstat_data   = ();
  my @Rstat_length = ();
  my @Rstat_names  = ();
  
  ## collect data for read.table string
  foreach my $sample (@{$self->data}){
    my ($basename,$dir,$ext) = fileparse($$sample{'bam'},qr/\.[^.]*/);
    
    push @Rstat_data, statsstring(@{$$sample{$whattodo}});
    push @Rstat_length, scalar(@{$$sample{$whattodo}});
    push @Rstat_names, "'$basename'",
    
  }

  my $data_string="summarydata<-list(stats=matrix(c(".join(",",@Rstat_data)."),5,".scalar(@Rstat_names)."), n=c(".join(",",@Rstat_length)."), names=c(".join(",",@Rstat_names)."))";

  ## produce box plot
  if(@Rstat_data){
    
    my $mn = "${whattodo}_stats.pdf";
    my $fn = file($self->outpath,$mn);
    
    $self->plot_bxplot($fn, $whattodo, $data_string);
  }
}

sub plot_bxplot{
  my ($self, $filename, $ylab, $datacommand_string) = @_;
  my ($bn,$odir,$ext) = fileparse($filename, qr /\..*/);
  my $rlibpath        = $self->rlib;
  #my $rlibpath        = '/usr/bin/R';
  $filename .= '.pdf' unless ($ext eq '.pdf');
  
  my $R = Statistics::R->new();
  $R->startR;
  $R->set('rlib', $rlibpath);
  $R->set('log_dir', $odir);
  $R->run("pdf('${filename}')") ;
  $R->run("$datacommand_string") ;
  $R->run("bxp(summarydata, medcol = 'red', ylab='$ylab', xlab='',las=3)") ;
  $R->run("dev.off()") ;
  $R->stopR;
}

sub statsstring{
  # usage: %h = %{stats(@a)};
  my @vals = sort {$a <=> $b} @_;
  my %stats = ();
  my @statstring = ();

  if(@vals){
    push @statstring, sprintf("%.2f", &min(\@vals));             ## min
    push @statstring, sprintf("%.2f", $vals[int(@vals/4)]);      ## 1.quartile
    if(@vals%2){
      push @statstring, $vals[int(@vals/2)];                     ## odd median
    }
    else{
      push @statstring, ($vals[int(@vals/2)-1] + $vals[int(@vals/2)])/2;  ## even median
    }
    push @statstring, sprintf("%.2f", $vals[int((@vals*3)/4)]);  ## 3.quartile
    push @statstring, sprintf("%.2f", &max(\@vals));             ## max
  }
  else{
    @statstring=qw/0 0 0 0 0/;
  }
  return(@statstring);
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

1;

__END__


=head1 NAME

Bio::ViennaNGS::BamStatSummary - Moose interface for analyzing,
summarizing and comparing BAM mapping statistics produced by
L<Bio::ViennaNGS::BamStat>

=head1 SYNOPSIS

  use Bio::ViennaNGS::BamStatSummary;

  $bamsummary->populate_data();
  $bamsummary->populate_countStat();
  $bamsummary->dump_countStat("csv");
  $bamsummary->make_BarPlot();

  $bamsummary->make_BoxPlot("data_edit" ) if( $bamsummary->control_edit   );
  $bamsummary->make_BoxPlot("data_clip" ) if( $bamsummary->control_clip   );
  $bamsummary->make_BoxPlot("data_match") if( $$bamsummary->control_match );
  $bamsummary->make_BoxPlot("data_qual" ) if( $bamsummary->control_qual   );


=head1 DESCRIPTION

This module provides a L<Moose> interface for processing the mapping
statistics of BAM files. Using the data structure populated by
L<Bio::ViennaNGS::BamStat>, it summarizes all data and compares
different BAM files. Output is generated both in CSV format and as
graphical representation of the results.  Internally, this modules
builds on L<Statistics::R>.


=head1 DEPENDENCIES

=over

=item L<Statistics::R>

=item L<Path::Class>

=item L<Tie::Hash::Indexed>

=item L<Moose>

=back

=head1 SEE ALSO

=over

=item L<Bio::ViennaNGS>

=item L<Bio::ViennaNGS::BamStat>

=back

=head1 AUTHORS

=over 

=item Fabian Amman E<lt>fabian@tbi.univie.ac.atE<gt>

=item Michael T. Wolfinger  E<lt>michael@wolfinger.euE<gt>

=back

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2015 by Michael T. Wolfinger

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut
