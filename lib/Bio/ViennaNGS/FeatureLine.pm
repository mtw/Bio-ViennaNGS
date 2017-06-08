# -*-CPerl-*-
# Last changed Time-stamp: <2017-06-08 18:34:41 michl>

package Bio::ViennaNGS::FeatureLine;

use version; our $VERSION = qv('0.17');
use namespace::autoclean;
use Moose;
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

