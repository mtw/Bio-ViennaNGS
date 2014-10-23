# -*-CPerl-*-
# Last changed Time-stamp: <2014-10-23 23:46:10 mtw>

package Bio::ViennaNGS::MinimalFeature;

use namespace::autoclean;

use Moose::Util::TypeConstraints;
subtype 'PlusOrMinus',
  as 'Str',
  where { /[\+\-\.]/ },
  message { "$_ is neither +/- or ."};
no Moose::Util::TypeConstraints;

use Moose;
has 'chr' => (
	      is  => 'rw',
	      isa => 'Str',
	      required => 1,
	     );

has 'start' => (
		is      => 'rw',
		isa     => 'Int',
		required => 1,
	       );

has 'end' => (
	      is      => 'rw',
	      isa     => 'Int',
	      required => 1,
	     );

has 'strand' => (
		 is      => 'rw',
		 isa     => 'PlusOrMinus',
		 default => '.',
		);

sub set_minimalFeature{
  my ($self,$chr,$start,$end,$strand) = @_;
  $self->chr($chr);
  $self->start($start);
  $self->end($end);
  $self->strand($strand);
}

no Moose;

# ============================================ #
package Bio::ViennaNGS::Feature;

use namespace::autoclean;
use Moose;
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

# ============================================ #
package Bio::ViennaNGS::FeatureChain;

use Tie::Hash::Indexed;
use namespace::autoclean;
use Moose;

has 'type' => (
	       is => 'rw',
	       isa => 'Str', # [mRNA,tRNA,rRNA.CDS,SJ,exon,intron,..]
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

