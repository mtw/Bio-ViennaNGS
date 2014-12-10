# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-10 15:53:45 mtw>

package Bio::ViennaNGS::FeatureChain;

use 5.12.0;
use version; our $VERSION = qv('0.11');

#use namespace::autoclean;
use Moose;
#use MooseX::InstanceTracking;

has 'type' => (
	       is => 'rw',
	       isa => 'Str', # [exon,intron,SJ,promoter,TSS,...]
	      );

has 'chain' => (
		is => 'rw',
		traits => ['Array'],
		isa => 'ArrayRef[Bio::ViennaNGS::Feature]',
		required => '1',
		builder => 'build_chain',
		auto_deref => 1,
		);

sub build_chain {
  my $self = shift;
  my $featurelist = shift; ## We expect features to be pre-sorted, So
                           ## I simply read in an arrayref of Feature
                           ## Objects
  return ($featurelist);
}

sub print_chain{
  my $self = shift;
  return 0 if (!$self->chain);
  my $out;
  foreach my $feature (@{$self->chain}){
    $out .= join("\t",
		 "chr".$feature->chromosome,
		 $feature->start,
		 $feature->end,
		 $feature->name,
		 $feature->score,
		 $feature->strand);
    $out .= "\t".$feature->extension if ($feature->extension);
    $out .= "\n";
  }
  return $out;
}

no Moose;

1;
