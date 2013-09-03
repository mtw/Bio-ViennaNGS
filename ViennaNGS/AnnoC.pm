package ViennaNGS::AnnoC;

use Exporter;
use strict;
use warnings;
use Data::Dumper;
use Bio::Tools::GFF;


our @ISA       = qw(Exporter);
our $VERSION   = '0.01';
our @EXPORT    = qw(parse_gff $feat $fstat);
our @EXPORT_OK = qw(%features %featstat);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

our %features  = ();
our %featstat  = ();
our $feat      = \%features;
our $fstat     = \%featstat;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

sub parse_gff {
  my ($i,$gffio,$feature,$gbkey);
  my $in_file = shift;

  $gffio = Bio::Tools::GFF->new(-file        => $in_file,
				-gff_version  => 3,
			       );
  $gffio->ignore_sequence(1);

  while($feature = $gffio->next_feature()) {
    my ($uid,$feat_name);
    my @name = my @id = my @gbkeys = ();

    next if ($feature->primary_tag() eq "exon");

    # 1. determine gbkey of the current feature
    @gbkeys = $feature->get_tag_values("gbkey");
    $gbkey  = $gbkeys[0];

    # 2. Get a unique ID for each feature
    if ($feature->has_tag('ID')){
      @id = $feature->get_tag_values('ID');
      $uid = $id[0]; # ID=id101
    }
    else {
      die "ERROR: feature '$gbkey' at pos. $feature->start does not have \'ID\' attribute\n";
    }

    # 3. assign parent's unique ID in case a parent record exists
    if ($feature->has_tag('Parent')){
      @id = $feature->get_tag_values('Parent');
      $uid = $id[0]; # ID=id101
    }

    # 4. Find a name for the current feature, use 'Name' or 'ID' attribute
    if ($feature->has_tag('Name')){
      @name = $feature->get_tag_values('Name');
      $feat_name = $name[0];
    }
    elsif ($feature->has_tag('ID')){
      @id = $feature->get_tag_values('ID');
      $feat_name = $id[0]; # ID=id101, use ID as feature name
    }
    else {
      die "ERROR: cannot set name for feature $feature->gbkey at pos. $feature->start\n";
    }

    unless (exists $features{$uid}) { # gene / ribosome_entry_site / etc.
      $features{$uid}->{start}     = $feature->start;
      $features{$uid}->{end}       = $feature->end;
      $features{$uid}->{strand}    = $feature->strand;
      $features{$uid}->{length}    = $feature->length;
      $features{$uid}->{seqid}     = $feature->seq_id;
      $features{$uid}->{score}     = $feature->score || 0;
      $features{$uid}->{gbkey}     = $gbkey;
      $features{$uid}->{name}      = $feat_name;
      $features{$uid}->{uid}       = $uid;
    }
    else { # CDS / tRNA / rRNA / etc
      #print "Overwriting \$f{$uid}, it was a $$f{$uid}->{gbkey} and will be a $gbkey\n";
      $features{$uid}->{gbkey} = $gbkey;  # gbkey for tRNA/ rRNA/ CDS etc
    }
  }

  # finally generate some statistics on features present in this annotation
  $featstat{total} = 0;
  foreach my $ft (keys %features){
    $featstat{total}++;
    my $key = $features{$ft}->{gbkey};
    unless (exists $featstat{$key}){
      $featstat{$key} = 0;
    }
    $featstat{$key} += 1;
  }
  $gffio->close();
 # print Dumper (\%features);
}


1;
__END__

=head1 NAME

ViennaNGS::AnnoC - Perl extension for converting annotation formats

=head1 SYNOPSIS

  use ViennaNGS::AnnoC;
  gff2bed($gff);

=head1 DESCRIPTION

gff2bed converts GFF3 to BED12. gff2bed expects the path to a GFF3
file and returns two hash references: $feat is a reference to a HOH
containing the raw annotation information for each feature found in
the GFF3 file. $fstat references a hash containing some statistics of
the features found in the GFF3 file.

=head2 EXPORT

Routines: gff2bed
Variables: feat fstat

=head1 AUTHOR

Michael T. Wolfinger, E<lt>michael@wolfinger.euE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 by Michael T. Wolfinger

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.12.4 or,
at your option, any later version of Perl 5 you may have available.


=cut
