# -*-CPerl-*-
# Last changed Time-stamp: <2015-02-24 13:56:09 mtw>

package Bio::ViennaNGS::ExtFeature;

use version; our $VERSION = qv('0.14');

use Moose;
with 'MooseX::Clone';
use MooseX::InstanceTracking;

extends 'Bio::ViennaNGS::Feature';

has 'extension' => (
		    is      => 'rw',
		    isa     => 'Str',
		    default => '',
		    predicate => 'has_extension',
		   );

sub set_ExtFeature{
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

1;
