# -*-CPerl-*-
# Last changed Time-stamp: <2015-10-22 13:01:02 mtw>
package Bio::ViennaNGS::FeatureIO;

use Moose;
use Carp;
use File::Slurp;
use Bio::ViennaNGS::Bed;
use Bio::ViennaNGS::Feature;
use Bio::ViennaNGS::BedGraphEntry;

#use Data::Dumper;

has 'file' => ( # file path
	       is => 'ro',
	       isa => 'Str',
	       predicate => 'has_file',
	       required => 1,
	      );

has 'filetype' => ( # BED6, BED12, GFF, GTF
		   is => 'ro',
		   isa => 'Str',
		   predicate => 'has_filetype',
		   required => 1,
		  );

has 'objecttype' => (
		 is => 'rw',
		 isa => 'Str', # BedGraph, ContainerFeature
		 predicate => 'has_objecttype',
		 required => 1,
		);

has 'data' => (
	       is => 'rw',
	       isa => 'ArrayRef',
	       default => sub { [] },
	       predicate => 'has_data',
	      );

sub BUILD { # call a parser method, depending on $self->objecttype
  my $self = shift;
  my $this_function = (caller(0))[3];

  #carp "INFO [$this_function]  now in BUILD";

  confess "ERROR [$this_function] \$self->file not available"
    unless ($self->has_file);
  confess "ERROR [$this_function] \$self->filetype not available"
    unless ($self->has_filetype);
  confess "ERROR [$this_function] \$self->objecttype not available"
    unless ($self->has_objecttype);

  if ($self->objecttype eq "BedGraph"){
    #carp "INFO  [$this_function] \$self->objecttype is BedGraph\n";
    $self->parse_bedgraph_file($self->file);
    return $self->data;
  }
  else{
    croak "ERROR [$this_function] Invalid type for \$self->objecttyp";
  }
}

sub parse_bedgraph_file{
  my ($self,$file) = @_;
  my $this_function = (caller(0))[3];
  my ($bg,$line,$entry,$chr,$start,$end,$val);

  $bg =  read_file( $file, array_ref => 1 ) ;
  foreach $line (@$bg){
    chomp($line);
    croak "ERROR [$this_function] cannot parse bedGraph input from $file"
      unless {$line =~ m/^([a-zA-Z0-9._]+)\t(\d+)\t(\d+)\t(-?\d+\.?\d*)$/};
    $chr = $1; $start = $2; $end = $3, $val = $4;
    #print "++ \$chr $chr ++ \$start $start ++ \$end $end ++ \$val $val\n";
    $entry = Bio::ViennaNGS::BedGraphEntry->new(chromosome => $chr,
						start => $start,
						end => $end,
						dataValue => $val,
					       );
    push @{$self->data}, $entry;
  }
}

no Moose;

1;
