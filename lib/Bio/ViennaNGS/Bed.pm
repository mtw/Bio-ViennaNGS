# -*-CPerl-*-
# Last changed Time-stamp: <2015-03-13 11:16:41 mtw>

package Bio::ViennaNGS::Bed;

use version; our $VERSION = qv('0.15_01');
use Carp;
use Moose;
use namespace::autoclean;

extends 'Bio::ViennaNGS::Feature';

has '+strand' => (
		  required => 1,
		 );

has 'thickStart' => (
		     is => 'rw',
		     isa => 'Int',
		     predicate => 'has_thickStart',
		    );

has 'thickEnd' => (
		   is => 'rw',
		   isa => 'Int',
		   predicate => 'has_thickEnd',
		  );

has 'itemRgb' => ( 
		  is => 'rw',
		  isa => 'Str',
		  lazy => 1,
		  default => '0',
		 );

has 'blockCount' => (
		     is => 'rw',
		     isa => 'Int',
		     predicate => 'has_blockCount',
		    );

has 'blockSizes' => (
		     is => 'rw',
		     isa => 'Value',
		    );

has 'blockStarts' => (
		      is => 'rw',
		      isa => 'Value',
		     );

has 'length' => (
		 is => 'rw',
		 isa => 'Int',
		 builder => '_build_length',
		 lazy => 1,
		 predicate => 'has_length',
		);


sub _build_length {
  my ($self) = @_;
  my ($i,$len) = (0)x2;
  my $this_function = (caller(0))[3];
  my @blockSizes = split (/,/ ,$self->blockSizes);
  my $bc = scalar @blockSizes;
  croak "ERROR [$this_function] invalid blockSount in BED12 line"
    unless ($bc == $self->blockCount);

  for ($i=0;$i<$bc;$i++){
    $len += $blockSizes[$i];
  }
  $self->length($len);
}

sub as_bed_line {
  my ($self,$n) = @_;
  my $this_function = (caller(0))[3];
  croak "ERROR [$this_function] argument of as_bed_line() must be 6 or 12"
    unless ( ($n == 6) | ($n == 12) );
  my $bed6= join ("\t",
		  $self->chromosome,
		  $self->start,
		  $self->end,
		  $self->name,
		  $self->score,
		  $self->strand,
		 );
  my $bed12 = join ("\t", $bed6,
		    $self->thickStart,
		    $self->thickEnd,
		    $self->itemRgb,
		    $self->blockCount,
		    $self->blockSizes,
		    $self->blockStarts,
		   );
  if ($n ==6){
    return $bed6;
  }
  else{
    return $bed12;
  }
}

# TODO:
# sub from_FeatureLine()
# sub to_FeatureLine()

__PACKAGE__->meta->make_immutable;

1;

__END__

=head1 NAME

Bio::ViennaNGS::Bed - Object-oriented interface for manipulation of
genomic interval data in BED format

=head1 SYNOPSIS

  use Bio::ViennaNGS::Bed;

  my $bedobject = Bio::ViennaNGS::Bed->new();

  # compute the length of a BED12 block
  $bedobject->_build_length();

  # dump an object as BED12 line
  $bedobject->as_bed_line(12);

=head1 DESCRIPTION

This module provides a Moose interface for storage and manipulation of
genomic interval data. It is primarily used as a convenience wrapper
for BED data with more generic L<Bio::ViennaNGS> classes for feature
annotation, such as L<Bio::ViennaNGS::MinimalFeature>,
L<Bio::ViennaNGS::Feature>, L<Bio::ViennaNGS::FeatureChain> and
L<Bio::ViennaNGS::FeatureLine>.

=head1 METHODS

=over

=item _build_length

Title   : _build_length

Usage   : C<$obj-E<gt>_build_length();>

Function: Compute the length of a BED12 interval block / line, i.e.
          the sum over the lengths of all intervals that make up a
          BED12 entry.

Args    :

Returns :

=item as_bed_line

Title   : as_bed_line

Usage   : C<$obj-E<gt>as_bed_line($bedtype);>

Function: Dump the contents of the object as BED6 or BED12 line.

Args    : C<$bedtype> can either be 6 or 12, determining BED6 or BED12
          output.

Returns : A (tab-separated) BED6 or BED12 line as string.

=back

=head1 DEPENDENCIES

=over

=item L<Moose>

=item L<Carp>

=item L<namespace::autoclean>

=back

=head1 AUTHORS

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2015 Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.


=cut
