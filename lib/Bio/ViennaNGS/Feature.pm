# -*-CPerl-*-
# Last changed Time-stamp: <2015-01-27 16:15:02 mtw>

package Bio::ViennaNGS::Feature;

use version; our $VERSION = qv('0.12_13');

#use namespace::autoclean;
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

has 'extension' => (
		    is      => 'rw',
		    isa     => 'Str',
		    default => '',
		   );

sub set_Feature{
  my ($self,$chromosome,$start,$end,$name,$score,$strand,$extension) = @_;
  $self->chromosome($chromosome);
  $self->start($start);
  $self->end($end);
  $self->name($name);
  $self->score($score);
  $self->strand($strand);
  $self->extension($extension);
}

no Moose;
#__PACKAGE__=>meta->make_immutable;

1;
