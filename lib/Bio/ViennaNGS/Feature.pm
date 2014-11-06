# -*-CPerl-*-
# Last changed Time-stamp: <2014-11-06 22:56:11 mtw>

package Bio::ViennaNGS::Feature;

#use namespace::autoclean;
use Moose;
use MooseX::InstanceTracking;
extends 'Bio::ViennaNGS::MinimalFeature';

has 'name' => (
	       is  => 'rw',
	       isa => 'Str',
	       required => '1',
	      );

has 'score' => (
		is => 'rw',
		isa => 'Value',
		default => '0',
	       );

has 'ext' => (
	      is      => 'rw',
	      isa     => 'Str',
	     );

no Moose;
#__PACKAGE__=>meta->make_immutable;

1;
