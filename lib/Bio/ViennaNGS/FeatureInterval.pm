# -*-CPerl-*-
# Last changed Time-stamp: <2015-10-21 11:44:21 mtw>

package Bio::ViennaNGS::FeatureInterval;

use version; our $VERSION = qv('0.16_01');
use namespace::autoclean;

use Moose;

has 'chromosome' => (
		     is  => 'rw',
		     isa => 'Str',
		     required => 1,
		     predicate => 'has_chromosome',
	     );

has 'start' => (
		is      => 'rw',
		isa     => 'Int',
		required => 1,
		predicate => 'has_start',
	       );

has 'end' => (
	      is      => 'rw',
	      isa     => 'Int',
	      required => 1,
	      predicate => 'has_end',
	     );


sub set_FeatureInterval{
  my ($self,$chr,$start,$end) = @_;
  $self->chromosome($chr);
  $self->start($start);
  $self->end($end);
}

no Moose;

1;

