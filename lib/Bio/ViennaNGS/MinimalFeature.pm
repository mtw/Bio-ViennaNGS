# -*-CPerl-*-
# Last changed Time-stamp: <2015-10-21 11:43:34 mtw>

package Bio::ViennaNGS::MinimalFeature;

use version; our $VERSION = qv('0.16_01');
use namespace::autoclean;

use Moose::Util::TypeConstraints;
subtype 'PlusOrMinus',
  as 'Str',
  where { /[\+\-\.]/ },
  message { "$_ is neither +/- nor ."};
no Moose::Util::TypeConstraints;

use Moose;
with 'MooseX::Clone';

extends 'Bio::ViennaNGS::FeatureInterval';

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

