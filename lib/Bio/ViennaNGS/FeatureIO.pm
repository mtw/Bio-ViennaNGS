# -*-CPerl-*-
# Last changed Time-stamp: <2017-06-10 19:07:22 michl>
package Bio::ViennaNGS::FeatureIO;

use Bio::ViennaNGS;
use Moose;
use Carp;
use File::Slurp;
use Bio::ViennaNGS::Bed;
use Bio::ViennaNGS::Feature;
use Bio::ViennaNGS::FeatureChain;
use Bio::ViennaNGS::BedGraphEntry;
use Data::Dumper;
use version; our $VERSION = version->declare("$Bio::ViennaNGS::VERSION");

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

has 'instanceOf' => (
		     is => 'rw',
		     isa => 'Str', # BedGraph, ContainerFeature
		     predicate => 'has_instanceOf',
		     required => 1,
		     writer => 'set_instanceOf',
		    );

has 'data' => (
	       is => 'rw',
	       isa => 'ArrayRef',
	       default => sub { [] },
	       traits => ['Array'],
	       predicate => 'has_data',
	       handles => {
			   all    => 'elements',
			   count  => 'count',
			   add    => 'push',
			   pop    => 'pop',
			  },
	      );

has '_entries' => ( # of elements in $self->data
		   is => 'rw',
		   isa => 'Int',
		   predicate => 'nr_entries',
		   init_arg => undef, # make this unsettable via constructor
		   builder => 'count_entries',
		   lazy => 1,
		  );

sub BUILD { # call a parser method, depending on $self->instanceOf
  my $self = shift;
  my $this_function = (caller(0))[3];
  my $type;

  confess "ERROR [$this_function] \$self->file not available"
    unless ($self->has_file);
  confess "ERROR [$this_function] \$self->filetype not available"
    unless ($self->has_filetype);
  confess "ERROR [$this_function] \$self->instanceOf not available"
    unless ($self->has_instanceOf);

  if ($self->filetype eq "BedGraph"){
    #carp "INFO  [$this_function] \$self->instanceOf is BedGraph\n";
    $self->parse_bedgraph_file($self->file);
    return $self->data;
  }
  elsif ($self->filetype =~ m/[Bb]ed6/){
    if($self->instanceOf eq "Feature"){
      #carp "INFO  [$this_function] \$self->instanceOf is Feature\n";
      $type=0; # ArrayRef of individual Feature objects
    }
    elsif ($self->instanceOf eq "FeatureChain"){
      #carp "INFO  [$this_function] \$self->instanceOf is FeatureChain\n";
      $type=1; # ArrayRef of FeatureChain objects, one per Feature object
    }
    elsif ($self->instanceOf eq "FeatureChainBlock"){
      #carp "INFO  [$this_function] \$self->instanceOf is FeatureChainBlock\n";
      $type=2; # ArrayRef of the entire block of Features (aka Bed12 from Bed6 block)
    }
    else{
      croak "ERROR [$this_function] Invalid type for \$self->instanceOf: $self->instanceOf";
    }
    $self->parse_bed6_file($self->file,$type);
    return $self->data;
  }
  elsif ($self->filetype =~ m/[Bb]ed12/){
    if($self->instanceOf eq "Bed"){
      #carp "INFO  [$this_function] \$self->instanceOf is Bed\n";
      $type=0; # ArrayRef of individual Bio::ViennaNGS::Bed objects
    }
    else {croak "ERROR [$this_function] currently only 'Bed' is a valid option for \$self->instance";}
    $self->parse_bed12_file($self->file,$type);
    return $self->data;
  }
  else{
    croak "ERROR [$this_function] Invalid type for \$self->filetyp: $self->filetype";
  }
  $self->count_entries();
}

sub count_entries {
  my $self = shift;
  my $cnt = scalar @{$self->data};
  $self->_entries($cnt);
}

sub parse_bedgraph_file{
  my ($self,$filename) = @_;
  my $this_function = (caller(0))[3];
  my ($file,$line,$entry,$chr,$start,$end,$val);

  $file =  read_file( $filename, array_ref => 1, chomp =>1 ) ;
  foreach $line (@$file){
    croak "ERROR [$this_function] cannot parse bedGraph input from $filename"
      unless {$line =~ m/^([a-zA-Z0-9._]+)\t(\d+)\t(\d+)\t(-?\d+\.?\d*)$/};
    $chr = $1; $start = $2; $end = $3, $val = $4;
    #print "++ \$chr $chr ++ \$start $start ++ \$end $end ++ \$val $val\n";
    $entry = Bio::ViennaNGS::BedGraphEntry->new(chromosome => $chr,
						start => $start,
						end => $end,
						dataValue => $val);
    push @{$self->data}, $entry;
  }
}

sub parse_bed6_file{
  my ($self,$file,$typ) = @_;
  my $this_function = (caller(0))[3];
  my ($line,$feat,$fc);
  $file = read_file( $file, array_ref => 1, chomp =>1 );

  if ($typ == 2){ # initialize an empty FeatureChain object
    $fc = Bio::ViennaNGS::FeatureChain->new(type => "feature");
   }

  #  print "********** in parse_bed6: typ= $typ ************\n";
  foreach $line (@$file){
    my @feat = split /\t/,$line;
    $feat = Bio::ViennaNGS::Feature->new(chromosome=>$feat[0],
					 start=>$feat[1],
					 end=>$feat[2],
					 name=>$feat[3],
					 score=>$feat[4],
					 strand=>$feat[5]);
    if($typ == 0){ # ArrayRef of individual Feature objects
      push @{$self->data}, $feat;
    }
    elsif ($typ == 1) { # ArrayRef of FeatureChain objects, one per Feature object
      $fc = Bio::ViennaNGS::FeatureChain->new(type => "feature",
					      chain => [$feat]);
      #      $fc->count_entries();
      push @{$self->data}, $fc;
    }
    elsif($typ == 2){
      $fc->add($feat);
      $fc->count_entries();
    }
    else{
      croak "ERROR [$this_function] don't know how to handle typ $typ";
    }
  } #end foreach
  if ($typ == 2) { push @{$self->data}, $fc; }
 # print Dumper($self);
}

sub parse_bed12_file{
  my ($self,$file,$typ) = @_;
  my $this_function = (caller(0))[3];
  my ($line,$feat,$fc);
  $file = read_file( $file, array_ref => 1, chomp =>1 );
  foreach $line (@$file){
    my @mcData = split /\t/,$line;
    if ($typ == 0){ # ArrayRef of Bio::ViennaNGS::Bed objects
      my $bo = Bio::ViennaNGS::Bed->new(chromosome   => $mcData[0],
					start        => $mcData[1],
					end          => $mcData[2],
					name         => $mcData[3],
					score        => $mcData[4],
					strand       => $mcData[5],
					thickStart   => $mcData[6],
					thickEnd     => $mcData[7],
					itemRgb      => $mcData[8],
					blockCount   => $mcData[9],
					blockSizes   => $mcData[10],
					blockStarts  => $mcData[11],
				       );
      push @{$self->data}, $bo;
    }
  }

}

no Moose;

1;

__END__

=head1 NAME

Bio::ViennaNGS::FeatureIO - Versatile I/O interface for Bio::ViennaNGS
feature annotation classes

=head1 SYNOPSIS

  use Bio::ViennaNGS::FeatureIO;

  # initialize a FeatureIO object from a Bed6 file
  my $data_bed = Bio::ViennaNGS::FeatureIO->new(
					       file => "file.bed6",
					       filetype => 'Bed6',
					       instanceOf => 'Feature',
					      );

  # initialize a FeatureIO object from a Bed12 file
  my $data_bed = Bio::ViennaNGS::FeatureIO->new(
					       file => "file.bed12",
					       filetype => 'Bed12',
					       instanceOf => 'Bed',
					      );

  # initialize a FeatureIO object from a bedGraph file
  my $obj = Bio::ViennaNGS::FeatureIO->new(file       => "file.bg",
                                           filetype   => "BedGraph",
                                           instanceOf => "BedGraph",
                                          );


=head1 DESCRIPTION

This module provides an object-oriented interface for easy
input/output operations on common feature annotation file formats. It
is - by design - a very generic module that stores all annotation data
within the C<$self-E<gt>data> ArrayRef. C<$self-E<gt>filetype>
specifies the file type to be processed. Currently parsing of I<Bed6>,
I<Bed12> and I<bedGraph> files is supported. C<$self-E<gt>instanceOf>
determines the object type of elements held in the C<$self-E<gt>data>
ArrayRef(s).

In case of parsing Bed6 data, C<$self-E<gt>instanceOf> can either be
C<Feature>, C<FeatureChain> or C<FeatureChainBlock>. While the first
causes C<$self-E<gt>data> to hold an ArrayRef to individual
L<Bio::ViennaNGS::Feature> objects, the second triggers creation of
individual L<Bio::ViennaNGS::FeatureChain> objects (each containing
exactly one feature interval, corresponding to individual Bed6
entries). A C<FeatureChainBlock> value to the C<$self-E<gt>instanceOf>
causes C<$self-E<gt>data> to hold an ArrayRef to a combined
L<Bio::ViennaNGS::FeatureChain> containing an entire block of
features. In the context of Bed annotation this corresponds to a
single Bed12 line (e.g. gene/transcript) that contains all individual
Bed6 features (e.g. exons). Evidently, this only makes sense if all
Bed6 features originate from the same chromosome and strand.

In case of parsing Bed12 data, currently only C<Bed> is supported for
C<$self-E<gt>instanceOf>, causing C<$self-E<gt>data> to hold an
ArrayRef to L<Bio::ViennaNGS::Bed> (aka Bed12) features. This will be
adjusted to L<Bio::ViennaNGS::FeatureChain> in the future.

In case of pasring bedGraph data, C<$self-E<gt>instanceOf> is ignored
and C<$self-E<gt>data> holds an ArrayRef to individual
L<Bio::ViennaNGS::BedGraph> objects.

=head1 METHODS

=over

=item parse_bedgraph_file

=item parse_bed6_file.

=item parse_bed12_file

These methods are used for object construction and should not be
called directly.

=back

=head1 DEPENDENCIES

=over

=item L<Bio::ViennaNGS::Bed>

=item L<Bio::ViennaNGS::Feature>

=item L<Bio::ViennaNGS::FeatureChain>

=item L<Bio::ViennaNGS::BedGraphEntry>

=item L<File::Slurp>

=item L<Carp>

=back

=head1 AUTHORS

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014-2016 Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.


=cut
