# -*-CPerl-*-
# Last changed Time-stamp: <2016-10-04 08:51:31 mtw>
package Bio::ViennaNGS::FeatureIO;

use Moose;
use Carp;
use File::Slurp;
use Bio::ViennaNGS::Bed;
use Bio::ViennaNGS::Feature;
use bio::ViennaNGS::FeatureChain;
use Bio::ViennaNGS::BedGraphEntry;
use Data::Dumper;

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
      carp "INFO  [$this_function] \$self->instanceOf is Feature\n";
      $type=0; # ArrayRef of individual Feature objects
    }
    elsif ($self->instanceOf eq "FeatureChain"){
      carp "INFO  [$this_function] \$self->instanceOf is FeatureChain\n";
      $type=1; # ArrayRef of FeatureChain objects, one per Feature object
    }
    elsif ($self->instanceOf eq "FeatureChainBlock"){
      carp "INFO  [$this_function] \$self->instanceOf is FeatureChainBlock\n";
      $type=2; # ArrayRef of the entire block of Features (aka Bed12 from Bed6 block)
    }
    else{
      croak "ERROR [$this_function] Invalid type for \$self->instanceOf: $self->instanceOf";
    }
    $self->parse_bed6_file($self->file,$type);
    return $self->data;
  }
  else{
    croak "ERROR [$this_function] Invalid type for \$self->filetyp: $self->filetype";
  }
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
      $fc->_count_entries();
      push @{$self->data}, $fc;
    }
    elsif($typ == 2){
      $fc->add($feat);
      $fc->_count_entries();
    }
    else{
      croak "ERROR [$this_function] don't know how to handle typ $typ";
    }
  } #end foreach
  if ($typ == 2) { push @{$self->data}, $fc; }
}

no Moose;

1;

__END__

=head1 NAME

Bio::ViennaNGS::FeatureIO - Versatile I/O interface for Bio::ViennaNGS
feature annotation classes

=head1 SYNOPSIS

  use Bio::ViennaNGS::FeatureIO;

  # create a new object and parse the contents of a bedGraph file
  my $obj = Bio::ViennaNGS::FeatureIO->new(file       => "file.bg",
                                           filetype   => "BedGraph",
                                           objecttype => "BedGraph",
                                          );



=head1 DESCRIPTION

This module provides an object-oriented interface for easy
input/output operations on common feature annotation file formats. It
is - by design - a very generic module that stores all annotation data
within the C<$self-E<gt>data> ArrayRef. C<$self-E<gt>objecttype>
determines the type of elements in C<$self-E<gt>data>.

Currently parsing of bedGraph files into an array of
L<Bio::ViennaNGS::BedGraphEntry> objects is supported.

=head1 METHODS

=over

=item parse_bedgraph_file

Title : parse_gff

Usage : C<$obj-E<gt>parse_bedgraph_file($bedgraph_file);>

Function: Parses bedGraph coverage data into C<$self-E<gt>data>.

Args : The full path to a bedGraph file

Returns : Nothing.

Notes : The bedGraph specification is available at
        L<http://genome.ucsc.edu/goldenpath/help/bedgraph.html>.

=back

=head1 DEPENDENCIES

=over

=item L<Bio::ViennaNGS::Bed>

=item L<Bio::ViennaNGS::Feature>

=item L<Bio::ViennaNGS::BedGraphEntry>

=item L<File::Slurp>

=item L<Carp>

=back

=head1 AUTHORS

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2014-2015 Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.0 or,
at your option, any later version of Perl 5 you may have available.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.


=cut
