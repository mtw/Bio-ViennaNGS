# -*-CPerl-*-
# Last changed Time-stamp: <2015-10-27 14:45:14 mtw>

package Bio::ViennaNGS::BedGraphEntry;

use version; our $VERSION = qv('0.16');

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
