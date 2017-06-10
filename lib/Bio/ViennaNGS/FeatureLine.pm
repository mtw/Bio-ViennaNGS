# -*-CPerl-*-
# Last changed Time-stamp: <2017-06-10 19:03:07 michl>

package Bio::ViennaNGS::FeatureLine;

use Bio::ViennaNGS;
use Moose;
use namespace::autoclean;
use version; our $VERSION = version->declare("$Bio::ViennaNGS::VERSION");

extends 'Bio::ViennaNGS::MinimalFeature';

has 'id' => (
	     is => 'rw',
	     isa => 'Str', # e.g. a transcript ID
	     required => '1',
	    );

has 'fc' => (
	     is => 'rw',
	     traits => ['Hash'],
	     isa => 'HashRef',
	     required => '1',
	     builder => '_build_fc',
	     auto_deref => '1',
	     handles => {
			 add => 'set',
			 elements => 'elements',
			 kv => 'kv',
			 count => 'count',
			 lookup => 'accessor',
			},
	    );



sub _build_fc {
  tie my %hash, 'Tie::Hash::Indexed';
  return \%hash; 
}

no Moose;

1;

