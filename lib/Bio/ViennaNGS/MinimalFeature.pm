# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-20 00:32:56 mtw>

package Bio::ViennaNGS::MinimalFeature;

use 5.12.0;
use version; our $VERSION = qv('0.12_07');
use namespace::autoclean;

use Moose::Util::TypeConstraints;
subtype 'PlusOrMinus',
  as 'Str',
  where { /[\+\-\.]/ },
  message { "$_ is neither +/- nor ."};
no Moose::Util::TypeConstraints;

use Moose;
with 'MooseX::Clone';

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

has 'strand' => (
		 is      => 'rw',
		 isa     => 'PlusOrMinus',
		 default => '.',
		 predicate => 'has_strand',
		);

sub set_minimalFeature{
  my ($self,$chr,$start,$end,$strand) = @_;
  $self->chromosome($chr);
  $self->start($start);
  $self->end($end);
  $self->strand($strand);
}

no Moose;

1;

