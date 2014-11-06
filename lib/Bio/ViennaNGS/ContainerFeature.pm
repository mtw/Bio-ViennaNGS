# -*-CPerl-*-
# Last changed Time-stamp: <2014-11-06 22:52:42 mtw>

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

1;
