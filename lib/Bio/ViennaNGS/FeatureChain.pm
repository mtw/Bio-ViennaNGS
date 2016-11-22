# -*-CPerl-*-
# Last changed Time-stamp: <2016-11-21 08:31:32 mtw>

package Bio::ViennaNGS::FeatureChain;

use version; our $VERSION = qv('0.17_01');
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
#  carp "INFO [$this_function] in 'before ";
  $self->sip( sub { $_[0]->start <=> $_[1]->start} )
};

#sub _get_start {
#  my $self = shift;
#  print ">>>>>>>>>>>>>>>>>>>>>\n";
#  my $this_function = (caller(0))[3];
#  carp "INFO [$this_function] in _get_start ";
#  print Dumper($self->chain);
#  my $element = {$self->chain}->[0];
#  confess "ERROR [$this_function] element is undef" unless (defined $element);
  #  print Dumper($element);
#  print "-----------------------\n";
#  return -99;
##  print Dumper($element);
##  $self->start($element->start);
#}

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
