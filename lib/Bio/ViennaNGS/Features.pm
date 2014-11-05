# -*-CPerl-*-
# Last changed Time-stamp: <2014-10-29 15:24:17 mtw>

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
	      predicate => 'has_chr',
	     );

has 'start' => (
		is      => 'rw',
		isa     => 'Int',
		required => 1,
		predicate => 'has_start',
	       );

has 'end' => (
	      is      => 'rw',
	      isa     => 'Int',
	      required => 1,
	      predicate => 'has_end',
	     );

has 'strand' => (
		 is      => 'rw',
		 isa     => 'PlusOrMinus',
		 default => '.',
		 predicate => 'has_strand',
		);

sub set_minimalFeature{
  my ($self,$chr,$start,$end,$strand) = @_;
  $self->chr($chr);
  $self->start($start);
  $self->end($end);
  $self->strand($strand);
}

#no Moose;

# ============================================ #
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

# ============================================ #
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

# ============================================ #
package Bio::ViennaNGS::FeatureLine;

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

# ============================================ #
package Bio::ViennaNGS::ContainerFeature;

use Moose;
extends 'Bio::ViennaNGS::MinimalFeature';
use Data::Dumper;

has '_uid' => ( # this one is populated in the BUILD method
	       is => 'rw',
	       isa => 'Str',
	       init_arg => undef,
	      );

has 'id' => (
	     is => 'rw',
	     isa => 'Str',
	     required => 1,
	     default => "unknown",
	    );

has 'type' => (
	       is => 'rw',
	       isa => 'Str',
	       required => 1,
	       default => 'generic',
	      );

has 'parent' => (
		 is => 'rw',
		 isa => 'Bio::ViennaNGS::ContainerFeature|Bio::ViennaNGS::FeatureIO',
		 required => 1,
		);

has 'line' => (
	       is => 'rw',
	       traits => ['Hash'],
	       isa => 'HashRef',
	       builder => '_build_line',
	       auto_deref => '1',
	       handles => {
			   add => 'set',
			   elements => 'elements',
			   kv => 'kv',
			   count => 'count',
			   lookup => 'accessor',
			  },
	      );

sub _build_line {
  tie my %hash, 'Tie::Hash::Indexed';
  return \%hash;
}

sub BUILD { # construct the 'uid' attribute
  my $self = shift;
  my $this_function = (caller(0))[3];

  confess "ERROR [$this_function] \$self->chr not available"
    unless ($self->has_chr);
  confess "ERROR [$this_function] \$self->start not available"
    unless ($self->has_start);
  confess "ERROR [$this_function] \$self->end not available"
    unless ($self->has_end);
  confess "ERROR [$this_function] \$self->strand not available"
    unless ($self->has_strand);
  $self->_uid(join('__',$self->chr,$self->start,$self->end,$self->strand));
}

no Moose;

# ============================================ #
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

