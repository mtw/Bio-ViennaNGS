# -*-CPerl-*-
# Last changed Time-stamp: <2017-06-08 19:13:07 michl>

package Bio::ViennaNGS::Peak;

use version; our $VERSION = qv('0.17');
use Moose;
use Carp;
use Data::Dumper;
use Path::Class;
use List::Util qw(sum sum0 min max first);
use Bio::ViennaNGS::Util qw(sortbed);

use namespace::autoclean;

has 'data' => (
	       is => 'ro',
	       isa => 'HashRef',
	       predicate => 'has_data',
	       default => sub { {} },
	      );

has 'region' => (
		 is => 'ro',
		 isa => 'HashRef',
		 predicate => 'has_region',
		 default => sub { {} },
		);

has 'peaks' => (
		 is => 'ro',
		 isa => 'HashRef',
		 predicate => 'has_peaks',
		 default => sub { {} },
	       );

has 'winsize' => (
		  is => 'ro',
		  isa => 'Int',
		  predicate => 'has_winsize',
		 );

has 'winsize' => (
		  is => 'ro',
		  isa => 'Int',
		  predicate => 'has_winsize',
		 );

has 'interval' => (
		   is => 'ro',
		   isa => 'Int',
		   predicate => 'has_interval',
		  );

has 'mincov' => (
		 is => 'ro',
		 isa => 'Int',
		 predicate => 'has_mincov',
		);

has 'length' => (
		 is => 'ro',
		 isa => 'Int',
		 predicate => 'has_length',
		);

has 'threshold' => (
		 is => 'ro',
		 isa => 'Value',
		 predicate => 'has_threshold',
		);

sub populate_data {
  my ($self,$filep,$filen) = @_;
  my $this_function = (caller(0))[3];
  my ($i,$element,$chr,$start,$end,$val,$lastend);
  my $have_lastend = 0;

  # parse [+] strand
  foreach $element (@{$filep->data}){
    unless (exists ${$self->region}{$element->chromosome}{pos}{start}){
      ${$self->region}{$element->chromosome}{pos}{start} = $element->start;
    }
    if($have_lastend == 1 && $element->start>$lastend){ # fill the gap with 0s
      for($i=$lastend;$i<$element->start;$i++){
	${$self->data}{$element->chromosome}{pos}[$i] = 0.;
      }
    }
    for ($i=$element->start;$i<$element->end;$i++){
      ${$self->data}{$element->chromosome}{pos}[$i]=$element->dataValue;
    }
    $lastend = $element->end;
    $have_lastend = 1;
  }

  $lastend = 0; # reset lastend before parsing [-] strand
  $have_lastend = 1; # needed to get the values parsed correctly below

  # parse [-] strand
  foreach $element (@{$filen->data}){
    if ($element->{dataValue} <= 0){
      $element->{dataValue} *= -1;
    }
    unless (exists ${$self->region}{$element->chromosome}{neg}{start}){
      ${$self->region}{$element->chromosome}{neg}{start} = $element->start;
    }
    if($have_lastend == 1 && $element->start>$lastend){ # fill the gap with 0s
      for($i=$lastend;$i<$element->start;$i++){
	${$self->data}{$element->chromosome}{neg}[$i] = 0.;
      }
    }
    for ($i=$element->start;$i<$element->end;$i++){
      ${$self->data}{$element->chromosome}{neg}[$i]=$element->dataValue;
    }
    $lastend = $element->end;
    $have_lastend = 1;
  }
  return;
}

sub raw_peaks {
  my ($self,$dest,$prefix,$log) = @_;
  my ($fn,$fn_u,$outfile,$maxidx,$have_pstart,$pstart,$pend,$chr);
  my ($winstart, $winend, $winsum, $mean, $lastmax, $from, $to);
  my ($index_of_maxwin);
  my $suffix = "rawpeaks.bed";
  my $this_function = (caller(0))[3];
  my $flex = 10;

  croak "ERROR [$this_function]: $dest does not exist\n"
    unless (-d $dest);

  if (defined $log){
    open(LOG, ">", $log) or croak "$!";
  }
  else {croak "ERROR [$this_function] \$log not defined";}

  $fn   = $prefix.".".$suffix;
  $fn_u = $prefix.".u.".$suffix;
  $outfile = file($dest,$fn_u);

  croak "ERROR [$this_function]: $self->data not available\n"
    unless ($self->has_data);

  croak "ERROR [$this_function]: $self->region not available\n"
    unless ($self->has_region);

  open(RAWPEAKS, ">", $fn_u) or croak $!;
  $lastmax = 0;
  foreach $chr (keys (%{$self->data})){

    # pos strand
    $maxidx      = $#{$self->data->{$chr}{pos}} ; # largest index in array
    $have_pstart = 0;
    $pstart      = -1;
    $pend        = -1;
    print LOG "processing chr $chr [+] strand $maxidx\n";
    croak "$chr has not been parsed in input data"
      unless (exists $self->region->{$chr}{pos}{start});

    for ($winstart=$self->region->{$chr}{pos}{start};$winstart<($maxidx-$self->winsize);$winstart+=$self->interval){
      # NOTE: does not handle the last window properly - > ignore this
      $winend = $winstart + $self->winsize - 1;
      $winsum = sum0 @{$self->data->{$chr}{pos}}[$winstart..$winend];
      $mean = $winsum/$self->winsize;
       #print LOG "**processing $winstart - $winend sum=$winsum mean=$mean lastmax=$lastmax\n";

      if ($mean > $lastmax){
	$lastmax = $mean;
	$index_of_maxwin = $winstart;
	unless ($have_pstart == 1) {
	  $pstart = $winstart;
	  $have_pstart = 1;
	  #print LOG "Peak start $chr:$winstart -- ";
	}
      }
      else {
	if ($mean > 0. && $mean <= $lastmax*$self->threshold ){ # peak stop criterion
	  $pend = $winend;
	  $have_pstart = 0;
	  $from = min($pstart,$pend);
	  $to   = max($pstart, $pend);
	  # now add +- $flex nt to raw peaks, required fro better peak boundary determination below
	  if ($from-$flex > 0){$from = $from - $flex;}
	  $to= $to+$flex; # we dont know the chromsize here, hence override it in case
	  #print LOG "$pend end\n";
	  my %pk = (
		    chr      => $chr,
		    start    => $from,
		    end      => $to,
		    strand   => 'pos',
		    maxindex => $index_of_maxwin,
		   );
	  push (@{ $self->peaks->{$chr} }, \%pk);
	  print LOG "**PEAK $chr [+] $from-$to max at $index_of_maxwin\n";
	  print RAWPEAKS "$chr\t$from\t$to\tRAW\t0\t+\n";

	  # reset parameters
	  $lastmax = 0;
	  $pstart = -1;
	  $pend   = -1;
	}
      }
    } # end for
    $lastmax = 0;

    # neg strand
    $maxidx      = $#{$self->data->{$chr}{neg}}; # largest index in array
    $have_pstart = 0;
    $pstart      = -1;
    $pend        = -1;
    print LOG "processing chr $chr [-] strand $maxidx\n";
    croak "$chr has not been parsed in input data"
      unless (exists $self->region->{$chr}{neg}{start});

    for ($winstart=$maxidx;$winstart>($self->region->{$chr}{neg}{start}+$self->winsize);$winstart-=$self->interval){
      # NOTE: does not handle the last window properly - > ignore this
      $winend = $winstart - $self->winsize + 1;
      $winsum = sum0 @{$self->data->{$chr}{neg}}[$winend..$winstart];
      $mean = $winsum/$self->winsize;
      #print "processing $winstart - $winend sum=$winsum mean=$mean lastmax=$lastmax\n";

      if ($mean > $lastmax){
	$lastmax = $mean;
	$index_of_maxwin = $winend;
	unless ($have_pstart == 1) {
	  $pstart = $winend;
	  $have_pstart = 1;
	  #print "Peak start $chr:$pstart -- ";
	}
      }
      else {
	if ($mean > 0. && $mean <= $lastmax*$self->threshold ){ # peak stop criterion
	  $pend = $winstart;
	  $have_pstart = 0;
	  $from = min($pstart,$pend);
	  $to   = max($pstart, $pend);
	  # now add +- $flex nt to raw peaks, required fro better peak boundary determination below
	  if ($from-$flex > 0){$from = $from - $flex;}
	  $to = $to+$flex; # we dont know the chromsize here, hence override it in case
	  #print "$pend end\n";
	  my %pk = (
		    chr    => $chr,
		    start  => $from,
		    end    => $to,
		    strand => 'neg',
		    maxindex => $index_of_maxwin,
		   );
	  push (@{ $self->peaks->{$chr} }, \%pk);
	  print LOG "**PEAK $chr [-] $from-$to max at $index_of_maxwin\n"; 
	  print RAWPEAKS "$chr\t$from\t$to\tRAW\t0\t-\n";

	  # reset parameters
	  $lastmax = 0;
	  $pstart = -1;
	  $pend   = -1;
	}
      }
    } # end for

    #print LOG Dumper($self->peaks);

  } # end foreach
  close(RAWPEAKS);
  close(LOG);
  sortbed($fn_u,$dest,$fn,1,$log);
}

sub final_peaks {
  my ($self,$dest,$prefix,$log) = @_;
  my ($fn,$fn_u,$outfile,$strand,$max,$idx,$val,$peak,$position,$pleft,$pright,$str);
  my $suffix = "candidatepeaks.bed";
  my $this_function = (caller(0))[3];

  croak "ERROR [$this_function]: $dest does not exist\n"
    unless (-d $dest);

  if (defined $log){
    open(LOG, ">>", $log) or croak "$!";
  }
  else {croak "ERROR [$this_function] \$log not defined";}

  $fn   = $prefix.".".$suffix;
  $fn_u = $prefix.".u.".$suffix;
  $outfile = file($dest,$fn_u);

  croak "ERROR [$this_function]: $self->data not available\n"
    unless ($self->has_data);

  croak "ERROR [$this_function]: $self->peaks not available\n"
    unless ($self->has_peaks);

  open(CANDIDATEPEAKS, ">", $fn_u) or croak $!;

  foreach my $chr (keys (%{$self->peaks})){
    print LOG "filtering candidate peaks in $chr\n";

    foreach $peak ( @{$self->peaks->{$chr}}){
      my $have_pleft = 0; # we have found a proper left end
      my $have_pright = 0; # we have found a proper right end
      $strand = $$peak{strand};
      print LOG ">>processing peak $$peak{start}-$$peak{end} [$strand] ";
      $max = max @{$self->data->{$chr}{$strand}}[$$peak{start}..$$peak{end}]; # max peak elevation
      print LOG "max is $max";
      if ($max < $self->mincov){print LOG " ..SKIPPING\n";next;}
      $idx = first { @{$self->data->{$chr}{$strand}}[$_] eq $max } $$peak{start}..$$peak{end}; # coordinate
      #print LOG Dumper(\$peak);
      print LOG " at $idx ..KEEPING\n";

      # process regions left/right of the maximum to identify peak boundaries
      # validity check first
      die "maximum at $idx is not within peak region"
	if ( $idx<$$peak{start} || $idx>$$peak{end} );

      # find left end
      for ($position = $idx; $position>$$peak{start}; $position--){
	$val =  ${$self->data->{$chr}{$strand}}[$position];
	print LOG "** VAL_L at pos $position $val - \$have_pleft=$have_pleft ";
	if ( $val < $max*$self->threshold){
	  $have_pleft = 1;
	  print LOG "exiting 'cause $val < ".$max*$self->threshold." - \$have_pleft=$have_pleft\n";
	  last;
	}
	print LOG "\n";
      }
      $pleft = $position;

      # find right end
      for ($position = $idx; $position<=$$peak{end}; $position++){
	$val =  ${$self->data->{$chr}{$strand}}[$position];
	print LOG "** VAL_R at pos $position $val - \$have_pright=$have_pright ";
	if ( $val < $max*$self->threshold){
	  $have_pright = 1;
	  print LOG "exiting 'cause $val < ".$max*$self->threshold."- \$have_pright=$have_pright\n";
	  last;
	}
	print LOG "\n";
      }
      $pright = $position;
      if ($strand eq "neg"){
	$str = '-';
      }
      else{
	$str = '+';
      }
      print LOG "# PEAK characterized at $chr:$pleft-$pright\n";

      if ($have_pleft == 1 && $have_pright == 1){ # print only if its a proper peak
	if ($pright-$pleft<=$self->length){ # filter peak length 
	  print CANDIDATEPEAKS "$chr\t$pleft\t$pright\tCANDIDATE\t$max\t$str\n";
	}
      }
    } # end foreach $peak
  } # end foreach $chr
  close(CANDIDATEPEAKS);
  sortbed($fn_u,$dest,$fn,1,$log);
}
__PACKAGE__->meta->make_immutable;

1;
__END__


=head1 NAME

Bio::ViennaNGS::Peak - An object oriented interface for characterizing
peaks in RNA-seq data

=head1 SYNOPSIS

  use Bio::ViennaNGS::Peak;

  # get an instance of Bio::ViennaNGS::peak
  my $peaks = Bio::ViennaNGS::Peak->new();

  # parse coverage for [+] and [-] strand from Bio::ViennaNGS::FeatureIO objects
  $peaks->populate_data($filep,$filen);

  # identify regions covered by RNA-seq signal ('raw peaks')
  $peaks->raw_peaks($dest,$prefix,$log);

  # characterize final peaks
  $peaks->final_peaks($dest,$prefix,$log);


=head1 DESCRIPTION

This module provides a L<Moose> interface for characterization of
peaks in RNA-seq coverage data.

=head1 METHODS

=over

=item populate_data

Title : populate_data

Usage : C<$obj-E<gt>populate_data($filep,$filen);>

Function : Parses RNA-seq coverage for positive and negative strand
            into C<@{$self-E<gt>data}>, a Hash of Arrays data structure.

Args : C<$filep> and C<$filen> are instances of L<Bio::ViennaNGS::FeatureIO>.

Returns : None.

Notes: The memory footprint of this method is rather high. It builds a
       Hash of Arrays data structure from L<Bio::ViennaNGS::FeatureIO>
       input objects of roughly the size of the underlying genome
       (chromosomes are hash keys, and there is an array containing
       coverage information for every genomic position referenced by
       hash values).

=item raw_peaks

Title : raw_peaks

Usage : C<$obj-E<gt>raw_peaks($dest,$prefix,$log);>

Function : This method identifies genomic regions ('raw peaks')
           covered by RNA-seq signal by means of a sliding window
           approach. RNA-seq coverage is read from
           C<@{$self-E<gt>data}> (which is populated by e.g. the
           C<populate_data> method). The sliding window approach
           processes [+] and [-] strand for all chromosomes in 5'
           -E<gt> 3' direction, whereby the mean value of each window
           is used as a representative for this window. Thereby both
           start and end coordinates, as well as position of the
           maximum elevation are identified. Here the end position of
           a covered region is defined as the coordinate of the window
           whose mean is less than a certain value
           (i.e. C<$self-E<gt>threshold> * peak maximum).

Raw peaks are stored in C<%{$self-E<gt>data}-E<gt>{peaks}>.

Args : C<$dest> contains the output path for results, C<$prefix> the
       prefix used for all output file names. C<$log> is the name of a
       log file, or undef if no logging is reuqired.

Returns : None. The output is a position-sorted BED6 file containing
          all raw peaks.

Notes : It is highly recommended to use I<normalized> input data in
        order to allow for multiple calls of this method with the same
        set of parameters on different samples.

=item final_peaks

Title : final_peaks

Usage : C<$obj-E<gt>final_peaks($dest,$prefix,$log);>

Function : This method characterizes final peaks from RNA-seq coverage
           found in C<%{$self-E<gt>data}-E<gt>{peaks}>. The latter is
           supposed to have been populated by C<$self-E<gt>raw_peaks>.

The procedure for finding final peaks is as follows: For each raw peak
found in C<%{$self-E<gt>data}-E<gt>{peaks}> the window of maximum
coverage is retrieved and a (second) sliding window approach is then
applied to regions both upstream and downstream of the maximum. Peak
boundaries are set at the position where the mean coverage of the
respective window is lower than C<$self-E<gt>threshold> * peak
maximum).

Peaks are reported if their total length (as determined by this
routine) is not longer than C<$self-E<gt>length>.

Args : C<$dest> contains the output path for results, C<$prefix> the
       prefix used for all output file names. C<$log> is the name of a
       log file, or undef if no logging is reuqired.

Returns : None. The output is a position-sorted BED6 file containing
          all candidate peaks.

Notes :

=back

=head1 DEPENDENCIES

=over

=item L<Moose>

=item L<Carp>

=item L<Path::Class>

=item L<List::Util>

=item L<namespace::autoclean>

=back

=head1 SEE ALSO

=over

=item L<Bio::ViennaNGS>

=item L<Bio::ViennaNGS::Util>

=back

=head1 AUTHOR

Michael T. Wolfinger, E<lt>michael@wolfinger.euE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2015-2017 by Michael T. Wolfinger

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

This software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

=cut
