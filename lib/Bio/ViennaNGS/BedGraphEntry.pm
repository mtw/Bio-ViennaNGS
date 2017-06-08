# -*-CPerl-*-
# Last changed Time-stamp: <2017-06-08 19:12:22 michl>

package Bio::ViennaNGS::BedGraphEntry;

use version; our $VERSION = qv('0.17');

use Moose;

extends 'Bio::ViennaNGS::FeatureInterval';

has 'dataValue' => (
		    is => 'rw',
		    isa => 'Num',
		    default => '0',
		    required => 1,
		    predicate => 'has_datavalue',
	       );

sub set_BedGraphEntry{
  my ($self,$chr,$start,$end,$dataValue) = @_;
  $self->chromosome($chr);
  $self->start($start);
  $self->end($end);
  $self->dataValue($dataValue);
}

no Moose;

1;
