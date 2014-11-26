# -*-CPerl-*-
# Last changed Time-stamp: <2014-11-10 19:34:35 fall>

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

has 'extension' => (
	      is      => 'rw',
	      isa     => 'Str',
	      default => '',
	     );

sub set_Feature{
  my ($self,$chromosome,$start,$end,$name,$score,$strand,$extension) = @_;
  $self->chromosome($chromosome);
  $self->start($start);
  $self->end($end);
  $self->name($name);
  $self->score($score);
  $self->strand($strand);
  $self->extension($extension);
}

no Moose;
#__PACKAGE__=>meta->make_immutable;

1;
