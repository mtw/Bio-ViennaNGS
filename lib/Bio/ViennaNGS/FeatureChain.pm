# -*-CPerl-*-
# Last changed Time-stamp: <2016-12-02 13:02:39 mtw>

package Bio::ViennaNGS::FeatureChain;

use version; our $VERSION = qv('0.17_02');
use Carp;
use Moose;
with 'MooseX::Clone';
use MooseX::InstanceTracking;
use Data::Dumper;

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
			    shift_chain => 'shift',
			    sip    => 'sort_in_place',
			   },
	       );

has 'start' => (
		is => 'ro',
		isa => 'Int',
		predicate => 'has_start',
		);

has '_entries' => (
		   is => 'rw',
		   isa => 'Int',
		   predicate => 'nr_entries',
		   init_arg => undef, # make this unsettable via constructor
		   builder => 'count_entries',
		   lazy => 1,
		  );

#has 'foo' => (
#	      is => 'ro',
#	      isa => 'Str',
#	      default => "AAA",
#	     );

sub BUILD {
  my $self = shift;
  my $this_function = (caller(0))[3];
 # carp"INFO [$this_function]";
  confess "ERROR [$this_function] \$self->chain not available"
    unless ($self->has_chain);
  $self->count_entries();
}

sub count_entries {
  my $self = shift;
  my $cnt = scalar @{$self->chain};
  $self->_entries($cnt);
}

before 'as_bed12_line' => sub {
  my $self = shift;
  my $this_function = (caller(0))[3];
  $self->sip( sub { $_[0]->start <=> $_[1]->start} )
};

sub sort_chain_ascending {
  my $self = shift;
  my $this_function = (caller(0))[3];
  $self->sip( sub { $_[0]->start <=> $_[1]->start} )
}

sub as_bed12_line{
  my ($self,$name,$score,$strand) = @_;
  my ($i,$chr,$start,$end,$feat,$bed12,$bsizes,$bstarts);
  my $count=0;
  my @blockSizes = ();
  my @blockStarts = ();
  # TODO check whether all features have the same chromosome id
  $chr   = @{$self->chain}[0]->chromosome;
  $start = @{$self->chain}[0]->start;
  $end   = @{$self->chain}[$#{$self->chain}]->end;
  unless (defined $name){$name=@{$self->chain}[0]->name;}
  unless (defined $score){$score=@{$self->chain}[0]->score;}
  unless (defined $strand){$strand=@{$self->chain}[0]->strand;}

  # TODO populate blockSizes and blockStarts
  for ($i=0;$i<=$#{$self->chain};$i++){
    $count++;
    $feat = @{$self->chain}[$i];
    push @blockSizes, eval($feat->end - $feat->start);
    push @blockStarts, ($feat->start - $start);;
  }
  $bsizes = join (",",@blockSizes);
  $bstarts = join (",", @blockStarts);
  $bed12 = join ("\t",$chr,$start,$end,$name,$score,$strand,$start,$end,"0",$count,$bsizes,$bstarts);
  return $bed12;
}

sub as_bed6_array{
  my $self = shift;
  $self->count_entries();
  my @bed6array=();
  return 0 unless ($self->has_chain);
  for (my $i=0;$i<$self->_entries;$i++){
    push @bed6array, join ("\t", 
			   @{$self->chain}[$i]->chromosome,
			   @{$self->chain}[$i]->start,
			   @{$self->chain}[$i]->end,
			   @{$self->chain}[$i]->name,
			   @{$self->chain}[$i]->score,
			   @{$self->chain}[$i]->strand);
  }
  return \@bed6array;
}

#sub clone {
#  my ( $self, %params ) = @_;
#  $self->meta->clone_object($self, %params);
#  return $self;
#}

no Moose;

1;

__END__

=head1 NAME

Bio::ViennaNGS::FeatureChain - Generic Moose wrapper class for
combined/linked genomic intervals, eg BED12 elements

=head1 SYNOPSIS

  use Bio::ViennaNGS::Feature;
  use Bio::ViennaNGS::FeatureChain;

  # get some new sequence features as instances of Bio::ViennaNGS::Feature
  my $start1 = 1100; my $start2 = 2345; my $start3 = 2987;
  my $end1 = 1346; my $end2 = 2544; my $end3 = 3076;
  my $name1 = "feat1"; my $name2 = "feat2"; my $name3 = "feat3";
  my $strand = "+";
  my $chr = "chr1";
  my $feat1 = Bio::ViennaNGS::Feature->(chromosome => $chr,
                                        start => $start1,
                                        end => $end1,
                                        name => $name1,
                                        strand => $strand,
                                        );
   my $feat2 = Bio::ViennaNGS::Feature->(chromosome => $chr,
                                        start => $start2,
                                        end => $end2,
                                        name => $name2,
                                        strand => $strand,
                                        );
   my $feat3 = Bio::ViennaNGS::Feature->(chromosome => $chr,
                                        start => $start3,
                                        end => $end3,
                                        name => $name3,
                                        strand => $strand,
                                        );

  # initialize a FeatureChain for two of these intervals
  my $fc = Bio::ViennaNGS::FeatureChain->new(type => "exon",
                                             chain => [$feat1,$feat3],
                                            );

  # append a genomic intervald to the chain
  $fc->add($feat2);

  # sort the chain in ascending order by start coordinates 
  $fc->sort_chain_ascending();

  # get the number of elements in the chain
  $fc->count_entries();


=head1 DESCRIPTION


=head1 METHODS

=over

=item sort_chain_ascending

=item as_bed12_line

=item as_bed6_array

=back

=head1 DEPENDENCIES

=over

=item L<Carp>

=back

=head1 AUTHORS

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014-2017 Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.


=cut
