# -*-CPerl-*-
# Last changed Time-stamp: <2014-11-06 22:48:49 mtw>
package Bio::ViennaNGS::FeatureIO;

use Moose;

has 'from' => ( # BED6, BED12, GFF, GTF
	       is => 'rw',
	       isa => 'Str',
	      );

has 'item' => (
	       is => 'rw',
	       isa => 'ArrayRef', # to Bio::ViennaNGS::ContainerFeature
	      );
no Moose;

1;
