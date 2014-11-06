# -*-CPerl-*-
# Last changed Time-stamp: <2014-11-06 22:54:36 mtw>

package Bio::ViennaNGS::FeatureChain;

use Tie::Hash::Indexed;
#use namespace::autoclean;
use Moose;
use MooseX::InstanceTracking;

has 'type' => (
	       is => 'rw',
	       isa => 'Str', # [exon,intron,SJ,promoter,TSS,...]
	      );

has 'chain' => (
		is => 'rw',
		traits => ['Hash'], # it's a HashRef to a tied hash
		isa => 'HashRef',
		required => '1',
		builder => '_build_chain',
		auto_deref => 1,
		handles => { #  Moose::Meta::Attribute::Native::Trait::Hash
			    add => 'set',
			    elements => 'elements',
			    kv => 'kv',
			    count => 'count',
			    lookup => 'accessor',
			   },
		);

sub _build_chain {
  tie my %hash, 'Tie::Hash::Indexed';
  return \%hash;
}

no Moose;

1;
