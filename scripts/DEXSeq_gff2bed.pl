#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-10-02 14:43:50 mtw>
#
# Convert DEXSeq flattened GFF file into BED (for easy UCSC browser
# integration)
#
# ***********************************************************************
# *  Copyright notice
# *
# *  Copyright 2013 Michael T. Wolfinger <michael@wolfinger.eu>
# *  All rights reserved
# *
# *  This program is free software: you can redistribute it and/or modify
# *  it under the terms of the GNU General Public License as published by
# *  the Free Software Foundation, either version 3 of the License, or
# *  (at your option) any later version.
# *
# *  This program is distributed in the hope that it will be useful,
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# *  GNU General Public License for more details.
# *
# *  You should have received a copy of the GNU General Public License
# *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
# *
# *  This copyright notice MUST APPEAR in all copies of the script!
# ***********************************************************************

use strict;
use warnings;
use Getopt::Long;
use File::Basename;
use Tie::File;
use Data::Dumper;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my $DEXSeq_gff = undef;
my $ext        = ".gff";
my $sortBed    = "sortBed";
my ($bednameuns,$bedname,$bednamer,$bn);
my %dex = ();
my @bedfile = ();

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

Getopt::Long::config('no_ignore_case');
&usage() unless GetOptions("i=s"             => \$DEXSeq_gff,
                           "-help"           => \&usage,
                           "v");

die "ERROR: no input GFF file provided" unless (defined $DEXSeq_gff);
die "ERROR: please use file extension .gff for input file" unless ($DEXSeq_gff =~ /\.gff$/);
$bn = basename($DEXSeq_gff, $ext);

# 1) parse DEXSeq flattened GFF file
open (DEXSEQGFF, "< $DEXSeq_gff") or die "Cannot open DEXSeq GFF file\n";

while(<DEXSEQGFF>){
  next if(/^##/); #ignore header
  chomp;
  my $key;
  my %attribs = ();
  my ($chr, $source, $type, $start, $end, $score,
      $strand, $phase, $attributes) = split("\t");
  #store nine columns in hash
  my %fields = (
		chr        => $chr,
		source     => $source,
		type       => $type,
		start      => $start,
		end        => $end,
		score      => $score,
		strand     => $strand,
		phase      => $phase,
		attributes => $attributes,
	       );

  # skip everysthing except lines with type=transcript
  next unless ($fields{type} eq 'exonic_part');
  my @add_attributes = split(";", $attributes);
  # store ids and additional information in second hash
  foreach my $attr ( @add_attributes ) {
    next unless $attr =~ /^\s*(.+)\s(.+)$/;
    my $c_type  = $1;
    my $c_value = $2;
    if($c_type  && $c_value){
      if(!exists($attribs{$c_type})){
	$attribs{$c_type} = [];
      }
      push(@{ $attribs{$c_type} }, $c_value);
    }
  }
  $fields{transcripts}        = $attribs{transcripts}->[0];
  $fields{exonic_part_number} = $attribs{exonic_part_number}->[0];
  $fields{gene_id}            = $attribs{gene_id}->[0];

  $fields{gene_id} =~ s/"//g;
  $fields{exonic_part_number} =~ s/"//g;
  $key = $fields{gene_id}."-E".$fields{exonic_part_number};

  if(!exists($dex{$key})){
    $dex{$key} = undef;
  }
   $dex{$key} = \%fields;
}
close(DEXSEQGFF);

$bednameuns = $bn.".uns.bed";
$bedname = $bn.".bed";
$bednamer = $bn."bedr";
open (BEDFILE, "> $bednameuns") or die "ERROR: cannot open $bednameuns for writing\n";
my $trackline = "track name=\"DEXSeq flat exons\" visibility=pack color=144,144,255\n";

foreach my $exon (keys %dex){
  my @bedline = ();
  $dex{$exon}->{start}--;
  @bedline= join ("\t", ($dex{$exon}->{chr},$dex{$exon}->{start},$dex{$exon}->{end},$exon,0,$dex{$exon}->{strand}));
  print BEDFILE "@bedline\n";
}

my $cmdl = $sortBed." -i ".$bednameuns." > ".$bedname;
system($cmdl);

tie @bedfile, 'Tie::File', $bedname;
unshift (@bedfile, $trackline);
untie @bedfile;

close(BEDFILE);

unlink($bednameuns);


#print Dumper(%dex);


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

sub usage {
 print <<EOF;

DEXSeq_GFF2BED.pl:  Convert DEXSeq flattened GFF to BED6

usage: $0 -i <DEXSEQ_FLAT_GFF> [options]
program specific options:                                   default:
 -i        <string> specify DEXSeq GFF file                 ($DEXSeq_gff)
 -help               print this information

EOF
exit;
}

