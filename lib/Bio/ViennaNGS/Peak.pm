# -*-CPerl-*-
# Last changed Time-stamp: <2015-10-16 15:52:02 mtw>

package Bio::ViennaNGS::Peak;

use version; our $VERSION = qv('0.16_01');
use Moose;
use Carp;
use Data::Dumper;
use Path::Class;
use File::Slurp;
use List::Util qw(sum sum0 min max first);
use Bio::ViennaNGS::Bed;
use Bio::ViennaNGS::Util qw(sortbed);
use Data::Dumper;

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



sub parse_coverage_bedgraph {
  my ($self,$filep,$filen) = @_;
  my $this_function = (caller(0))[3];
  my ($i, $bg_pos, $bg_neg, $chr, $start, $end, $val, $lastend, $line);
  my $have_lastend = 0;

  croak "ERROR [$this_function] bedGraph input file [+]  $filep not available\n"
    unless (-e $filep);
  croak "ERROR [$this_function] bedGraph input file [-]  $filen not available\n"
    unless (-e $filen);

  $bg_pos = read_file( $filep, array_ref => 1 ) ;
  $bg_neg = read_file( $filen, array_ref => 1 ) ;

  # parse [+] strand
  foreach $line (@$bg_pos) {
    chomp($line);
    croak "ERROR [$this_function] cannot parse bedGraph [+] input"
      unless {$line =~ m/^([a-zA-Z0-9._]+)\t(\d+)\t(\d+)\t(\d+\.?\d*)$/};
    $chr = $1; $start = $2; $end = $3, $val = $4;
    unless (exists ${$self->region}{$chr}{pos}{start}){
      ${$self->region}{$chr}{pos}{start} = $start;
    }
    if($have_lastend == 1 && $start>$lastend){ # fill the gap with 0s
      for($i=$lastend;$i<$start;$i++){
	${$self->data}{$chr}{pos}[$i] = 0.;
      }
    }
    for ($i=$start;$i<$end;$i++){
      ${$self->data}{$chr}{pos}[$i]=$val;
    }
    $lastend = $end;
    $have_lastend = 1;
  }

  $lastend = 0; # reset lastend before parsing [-] strand
  $have_lastend = 1; # needed to get the values parsed correctly below

  # parse [-] strand
  foreach $line (@$bg_neg) {
    chomp($line);
    die "ERROR [$this_function] cannot parse bedGraph [-] input"
      unless {$line =~ m/^([a-zA-Z0-9._]+)\t(\d+)\t(\d+)\t(-?\d+\.?\d*)$/};
    $chr = $1; $start = $2; $end = $3, $val = $4;
    if ($val <= 0){
      $val *= -1;
    }
    unless (exists ${$self->region}{$chr}{neg}{start}){
      ${$self->region}{$chr}{neg}{start} = $start;
    }
    if($have_lastend == 1 && $start>$lastend){ # fill the gap with 0s
      for($i=$lastend;$i<$start;$i++){
	${$self->data}{$chr}{neg}[$i] = 0.;
      }
    }
    for ($i=$start;$i<$end;$i++){
      ${$self->data}{$chr}{neg}[$i]=$val;
    }
    $lastend = $end;
    $have_lastend = 1;
  }
}

sub raw_peaks {
  my ($self,$dest,$prefix,$log) = @_;
  my ($fn,$fn_u,$outfile,$maxidx,$have_pstart,$pstart,$pend,$chr);
  my ($winstart, $winend, $winsum, $mean, $lastmax, $from, $to);
  my ($index_of_maxwin);
  my $suffix = "rawpeaks.bed";
  my $this_function = (caller(0))[3];

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
      # print "processing $winstart - $winend sum=$winsum mean=$mean lastmax=$lastmax\n";

      if ($mean > $lastmax){
	$lastmax = $mean;
	$index_of_maxwin = $winstart;
	unless ($have_pstart == 1) {
	  $pstart = $winstart;
	  $have_pstart = 1;
	  #print "Peak start $chr:$winstart -- ";
	}
      }
      else {
	if ($mean > 0. && $mean <= $lastmax*$self->threshold ){ # peak stop criterion
	  $pend = $winend;
	  $have_pstart = 0;
	  $from = min($pstart,$pend);
	  $to   = max($pstart, $pend);
	  #print "$pend end\n";
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
      $max = max @{$self->data->{$chr}{$strand}}[$$peak{start}..$$peak{end}]; # max peak elevation
      next if ($max < $self->mincov);
      $idx = first { @{$self->data->{$chr}{$strand}}[$_] eq $max } $$peak{start}..$$peak{end}; # coordinate
      #print LOG Dumper(\$peak);
      print LOG "max in window $$peak{start}..$$peak{end} is $max at $idx\n";

      # process regions left/right of the maximum to identify peak boundaries
      # validity check first
      die "maximum at $idx is not within peak region"
	if ( $idx<$$peak{start} || $idx>$$peak{end} );

      # find left end
      for ($position = $idx; $position>$$peak{start}; $position--){
	$val =  ${$self->data->{$chr}{$strand}}[$position];
	#print "* VAL_L at pos $position $val - \$have_pleft=$have_pleft ";
	if ( $val < $max*$self->threshold){
	  $have_pleft = 1;
	  #print "exiting 'cause $val < ".$max*$threshold." - \$have_pleft=$have_pleft\n";
	  last;
	}
	#print "\n";
      }
      $pleft = $position;

      # find right end
      for ($position = $idx; $position<=$$peak{end}; $position++){
	$val =  ${$self->data->{$chr}{$strand}}[$position];
	#print "* VAL_R at pos $position $val - \$have_pright=$have_pright ";
	if ( $val < $max*$self->threshold){
	  $have_pright = 1;
	  #print "exiting 'cause $val < ".$max*$threshold."- \$have_pright=$have_pright\n";
	  last;
	}
	#print "\n";
      }
      $pright = $position;
      if ($strand eq "neg"){
	$str = '-';
      }
      else{
	$str = '+';
      }
      #print LOG "# PEAK characterized at $chr:$pleft-$pright\n";

      if ($have_pleft == 1 && $have_pright == 1){ # print only if its a proper peak
	if ($pright-$pleft<=$self->length){ # filter peak length 
	  print CANDIDATEPEAKS "$chr\t$pleft\t$pright\tCANDIDATE\t$max\t$str\n";
	}
      }
    } # end foreach $peak
  } # end foreach $chr
  close(CANDIDATEPEAKS);
}
__PACKAGE__->meta->make_immutable;

1;
__END__


=head1 NAME

Bio::ViennaNGS::Peak - An object oriented interface for characterizing
peaks in RNA-seq data


=head1 SYNOPSIS

  use Bio::ViennaNGS::Peak;

  my $peaks = Bio::ViennaNGS::Peak->new();

  # parse coverage from a BedGraph file
  $peaks->parse_coverage_bedgraph($filep,$filen);

  # identify regions covered by RNA-seq signal ('raw peaks')
  $peaks->raw_peaks();

  # characterize final peaks
  $peaks->final_peaks();


=head1 DESCRIPTION

This module provides a L<Moose> interface for characterization of
peaks in RNA-seq signal, provided via bedGraph input.

=head1 METHODS

=over

=item parse_coverage_bedgraph

Title : parse_coverage_bedgraph

Usage : C<$obj-E<gt>parse_coverage_bedgraph($filep,$filen);>

Function : Parses RNA-seq coverage from bedGraph input files for
           positive ($filep) and negative ($filen) strand into a Hash
           of Arrays data structure (C<@{$self-E<gt>data}>).

Args : C<$file> is a bedGraph input file.

Returns :

Notes :


=back

=head1 DEPENDENCIES

=over

=item L<Moose>

=item L<Carp>

=item L<Path::Class>

=item L<File::Slurp>

=item L<List::Util>

=item L<namespace::autoclean>

=back

=head1 SEE ALSO

=over
