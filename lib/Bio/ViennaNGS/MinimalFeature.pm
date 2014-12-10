# -*-CPerl-*-
# Last changed Time-stamp: <2014-11-10 19:40:06 fall>

package Bio::ViennaNGS::MinimalFeature;

use namespace::autoclean;

use Moose::Util::TypeConstraints;
subtype 'PlusOrMinus',
  as 'Str',
  where { /[\+\-\.]/ },
  message { "$_ is neither +/- nor ."};
no Moose::Util::TypeConstraints;

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

