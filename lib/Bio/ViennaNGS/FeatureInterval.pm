# -*-CPerl-*-
# Last changed Time-stamp: <2017-06-08 19:15:52 michl>

package Bio::ViennaNGS::FeatureInterval;

use version; our $VERSION = qv('0.17');
use namespace::autoclean;
use Carp;
use Moose;

has 'chromosome' => (
		     is  => 'ro',
		     isa => 'Str',
		     required => 1,
		     predicate => 'has_chromosome',
	     );

has 'start' => (
		is      => 'ro',
		isa     => 'Int',
		required => 1,
		predicate => 'has_start',
	       );

has 'end' => (
	      is      => 'ro',
	      isa     => 'Int',
	      required => 1,
	      predicate => 'has_end',
	     );

has '_length' => ( # length of interval
		   is => 'rw',
		   isa => 'Int',
		   predicate => 'length',
		   init_arg => undef, # make this unsettable via constructor
		 );

sub BUILD { # call a parser method, depending on $self->instanceOf
  my $self = shift;
  my $this_function = (caller(0))[3];

  confess "ERROR [$this_function] \$self->end must be >= than \$self->start"
    unless ($self->end >= $self->start);
  my $len = $self->end-$self->start;
  $self->_length($len);
}

no Moose;

1;

__END__

=head1 NAME

Bio::ViennaNGS::FeatureInterval - A Moose class for unstranded, elementary
genomic intervals.

=head1 SYNOPSIS

  use Bio::ViennaNGS::FeatureInterval;

  my $obj = Bio::ViennaNGS::FeatureInterval->new(chromosome => "chr1",
                                                 start => "1200",
                                                 end => "4300",
                                                );

  my $len = $obj->_length();

=head1 DESCRIPTION

This module provides an object-oriented interface for storing
unstranded, elementary genomic intervals characterized via chromosome,
start and end position. As such, it can be regarded a simple wrapper
for BED3 elements.

This is the base class for L<Bio::ViennaNGS::MinimalFeature>.

=head1 SEE ALSO

=over

=item L<Bio::ViennaNGS>

=item L<Bio::ViennaNGS::MinimalFeature>

=back

=head1 AUTHORS

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014-2017 Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

=cut
