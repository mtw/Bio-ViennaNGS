# -*-CPerl-*-
# Last changed Time-stamp: <2016-10-03 16:41:12 mtw>

package Bio::ViennaNGS::FeatureChain;

use version; our $VERSION = qv('0.17_01');
use Carp;
use Moose;
with 'MooseX::Clone';
use MooseX::InstanceTracking;

has 'type' => (
	       is => 'rw',
	       isa => 'Str', # [exon,intron,SJ,promoter,TSS,...]
	      );

has 'chain' => (
		is => 'rw',
		traits => ['Array', 'Clone' => {to=>'ArrayRef'}],
		isa => 'ArrayRef[Bio::ViennaNGS::Feature]',
		default => sub { [] },
		predicate => 'has_chain',
		auto_deref => 1,
		handles => {
			    all    => 'elements',
			    count  => 'count',
			    add    => 'push',
			    pop    => 'pop',
			   },
	       );

has '_entries' => (
		   is => 'rw',
		   isa => 'Int',
		   predicate => 'nr_entries',
		   init_arg => undef, # make this unsettable via constructor
		  );

sub BUILD {
  my $self = shift;
  my $this_function = (caller(0))[3];
  confess "ERROR [$this_function] \$self->chain not available"
    unless ($self->has_chain);
  $self->_entries( scalar @{$self->chain});
}

sub print_chain{
  my $self = shift;
  return 0 unless ($self->has_chain);
  my $out;
  foreach my $feature (@{$self->chain}){
    $out .= join("\t",
		 "chr".$feature->chromosome,
		 $feature->start,
		 $feature->end,
		 $feature->name,
		 $feature->score,
		 $feature->strand);
    $out .= "\n";
  }
  return $out;
}
  
  sub as_bed12_line{
    return;
  }
  
  #sub clone {
  #  my ( $self, %params ) = @_;
  #  $self->meta->clone_object($self, %params);
  #  return $self;
  #}
  
no Moose;

1;
