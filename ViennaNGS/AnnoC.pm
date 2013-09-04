#!/usr/bin/perl -w
# -*-CPerl-*-
# Last changed Time-stamp: <2013-09-04 15:33:34 mtw>
#
#
# ***********************************************************************
# *  Copyright notice
# *
# *  Copyright 2013 Michael Thomas Wolfinger <michael@wolfinger.eu>
# *  All rights reserved
# *
# * This library is free software; you can redistribute it and/or modify
# * it under the same terms as Perl itself, either Perl version 5.12.4 or,
# * at your option, any later version of Perl 5 you may have available.
# *
# *  This program is distributed in the hope that it will be useful,
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# *
# ***********************************************************************

package ViennaNGS::AnnoC;

use Exporter;
use strict;
use warnings;
use Data::Dumper;
use Bio::Tools::GFF;

our @ISA       = qw(Exporter);
our $VERSION   = '0.02';
our @EXPORT    = qw(parse_gff feature_summary $feat $fstat);
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
  my ($i,$in_file,$gffio,$feature,$gbkey);
  $in_file = shift;

  $gffio = Bio::Tools::GFF->new(-file        => $in_file,
				-gff_version  => 3,
			       );
  $gffio->ignore_sequence(1);
  if (my $header = $gffio->next_segment() ){
    $featstat{accession}= $header->display_id();
  }
  else{
    warn "INFO parse_gff: could not parse GFF header\n";
  }

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
      $features{$uid}->{gbkey} = $gbkey;  # gbkey for tRNA/ rRNA/ CDS etc
    }
  }

  # finally generate some statistics on features present in this annotation
  $featstat{total} = 0;
  $featstat{origin} = "ViennaNGS::AnnoC::parse_gff version ".$VERSION;
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

# feature_summary($fstat,$path)
# Print summary of $featstat hash
#
# ARG1: reference to $featstat hash
# ARG2: path for output file summary.txt
sub feature_summary {
  my ($fn, $summary_name);
  my ($summary, $wd) = @_;
  #print Dumper($summary);
  $fn = $$summary{accession}.".summary.txt";
  $wd .= "/" unless ($wd =~ /\/$/);
  $summary_name = $wd.$fn;
  open my $smry, ">", $summary_name or die $!;

  print $smry "Accession\t $$summary{accession}\n";
  print $smry "Origin   \t $$summary{origin}\n";
  foreach my $ft (sort keys %$summary){
    next if ($ft =~ /total/ || $ft =~ /accession/ || $ft =~ /origin/);
    print $smry "$ft\t$$summary{$ft}\n";
  }
  print $smry "Total\t$$summary{total}\n";
  close $smry;
}

1;
__END__

=head1 NAME

ViennaNGS::AnnoC - Perl extension for converting sequence annotation formats

=head1 SYNOPSIS

  use ViennaNGS::AnnoC;

  parse_gff($gff3_file);
  feature_summary($fstat,$path)

=head1 DESCRIPTION

parse_gff parses GFF3 annotation files. The GFF3 specification is
available at http://www.sequenceontology.org/resources/gff3.html
parse_gff expects the path to a GFF3 file as argument and returns two
hash references: $feat is a reference to a HOH containing the raw
annotation information for each feature found in the GFF3 file. $fstat
references a hash containing summary statistics of the features found
in the GFF3 file.

feature_summary generates a summary file for all features parsed by
parse_gff. feature_summary expects two arguments, first a refence to
the summary hash generated by parse_gff and second the path where the
summary.txt output should be written.

=head2 EXPORT

Routines: parse_gff feature_summary
Variables: feat fstat

=head1 AUTHOR

Michael Thomas Wolfinger, E<lt>michael@wolfinger.euE<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 by Michael Thomas Wolfinger

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.12.4 or,
at your option, any later version of Perl 5 you may have available.


=cut
