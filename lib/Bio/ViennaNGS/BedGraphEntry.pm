# -*-CPerl-*-
# Last changed Time-stamp: <2017-06-10 18:21:46 michl>

package Bio::ViennaNGS::BedGraphEntry;

use Bio::ViennaNGS;
use Moose;
use version; our $VERSION = version->declare("$Bio::ViennaNGS::VERSION");


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
