# -*-CPerl-*-
# Last changed Time-stamp: <2014-11-06 22:57:05 mtw>

package Bio::ViennaNGS::MinimalFeature;

use namespace::autoclean;

use Moose::Util::TypeConstraints;
subtype 'PlusOrMinus',
  as 'Str',
  where { /[\+\-\.]/ },
  message { "$_ is neither +/- or ."};
no Moose::Util::TypeConstraints;

use Moose;
has 'chr' => (
	      is  => 'rw',
	      isa => 'Str',
	      required => 1,
	      predicate => 'has_chr',
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
  $self->chr($chr);
  $self->start($start);
  $self->end($end);
  $self->strand($strand);
}

#no Moose;

1;

