# -*-CPerl-*-
# Last changed Time-stamp: <2015-02-24 13:56:31 mtw>

package Bio::ViennaNGS::Feature;

use version; our $VERSION = qv('0.14');

use Moose;
with 'MooseX::Clone';
use MooseX::InstanceTracking;

extends 'Bio::ViennaNGS::MinimalFeature';

has 'name' => (
	       is  => 'rw',
	       isa => 'Str',
	       required => '1',
	       predicate => 'has_name',
	      );

has 'score' => (
		is => 'rw',
		isa => 'Value',
		default => '0',
		predicate => 'has_score',
	       );

sub set_Feature{
  my ($self,$chromosome,$start,$end,$name,$score,$strand) = @_;
  $self->chromosome($chromosome);
  $self->start($start);
  $self->end($end);
  $self->name($name);
  $self->score($score);
  $self->strand($strand);
}

no Moose;

1;
