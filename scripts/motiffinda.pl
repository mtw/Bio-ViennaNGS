#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-10-09 17:08:31 mtw>
#
# Find motifs in annotated sequences. The motif is provided as regex
# via the command line
#
# usage: motiffinda.pl -motif <REXGEX> -gff <GFFFILE> -fa <FASTAFILE>
#                      -offset <OFFSET> -inframe
#
# The motif must be privided within braces, because we use $1 to obtain
# the position of the motif within the query sequence.
#
# ./motiffinda.pl -motif "(ACATG\w{4,13}ACATG)" -gff NC_000913.gff
#                 -fa NC_000913.fna -offset 1
#
# ***********************************************************************
# *  Copyright notice
# *
# *  Copyright 2014 Michael T. Wolfinger <michael@wolfinger.eu>
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
#

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Bio::ViennaNGS::Fasta;
use Bio::ViennaNGS::AnnoC qw(&parse_gff &get_fasta_ids @fastaids $feat $fastadb);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my ($gff_in,$fa_in,$obj,$fastaO);
my $motif      = undef;
my $inframe    = 0;
my $offset     = 0;
my @fastaIDs   = ();

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

Getopt::Long::config('no_ignore_case');
&usage() unless GetOptions("gff=s"      => \$gff_in,
                           "fa=s"       => \$fa_in,
                           "motif=s"    => \$motif,
			   "inframe"    => \$inframe,
			   "offset=i"   => \$offset,
                           "-help"      => \&usage,
                           "v");
unless ($gff_in =~ /^\//) {$gff_in = "./".$gff_in;}
die "$gff_in not found\n" unless (-f $gff_in);

unless ($fa_in =~ /^\//) {$fa_in = "./".$fa_in;}
die "$fa_in not found\n" unless (-f $fa_in);

unless (defined $motif) {die "please provide motif as regex\n";}

$fastaO = Bio::ViennaNGS::Fasta->new(fa=>$fa_in);

parse_gff($gff_in);
find_motifs($fastaO);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

# loop over all features parsed from GFF3 and see if we have the motif
sub find_motifs {
  my ($fo) = @_;
  my ($id,$start,$stop,$strand,$name,$seq,$length,$chr);

  foreach $id (keys %$feat){
    my $have_motif = 0;
    my @fm = ();  # AoH for storing found motifs and their position

    # TODO let user choose where to look for the mtoif 
    next unless ($$feat{$id}->{gbkey} eq "CDS");
    $start  = $$feat{$id}->{start};
    $stop   = $$feat{$id}->{end};
    $strand = $$feat{$id}->{strand};
    $name   = $$feat{$id}->{name};
    $length = $$feat{$id}->{length};
    $chr    = $$feat{$id}->{seqid};

    # check if we can provide sequence information for this $chr
    unless ($fo->has_sequid($chr)){
      my $_ids= $fo->fastaids;
      print STDERR "\'$chr\' is not a valid Fasta ID. ";
      print STDERR "The provided Fasta file provides the following IDs:\n";
      print STDERR join " | ", @$_ids;
      print STDERR "\n";
      die;
    }

    $seq    = $fo->stranded_subsequence($chr,$start,$stop,$strand);
    #print Dumper($$feat{$id});
    #print Dumper(\$seq);

    #  $seq    =~ s/T/U/g; # make RNA from DNA
    while ($seq =~ m/$motif/g) {
      my $p = pos($seq)-length($1)+1;
      $have_motif++;
      push (@fm, {motif=>$1,
		  pos=>$p,
		  name=>$name,
		  strand=>$strand,
		  start=>$start,
		  stop=>$stop},
	   );
    }
    
    print Dumper(\@fm);
    for (my $i=0;$i<$have_motif;$i++) {
      next unless (($inframe == 1) && (($fm[$i]->{pos}+$offset)%3==0));
      # ATG in frame
      print ">$fm[$i]->{name}|$fm[$i]->{strand}|rel $fm[$i]->{pos}|";
      if ($fm[$i]->{strand} == 1){
	print "abs ".eval($fm[$i]->{pos}+$fm[$i]->{start}-1)."\n";
      }
      else {
	print "abs ".eval($fm[$i]->{stop}-$fm[$i]->{pos}+1)."\n";
      }
      print "$fm[$i]->{motif}\n";
    }
  }
}


sub usage {
  print <<EOF;
Find  motifs in annotated sequences

usage: $0 [options]
program specific options:                                    default:
 -motif    <string>  motif to search for (as regex)           ($motif)
 -gff      <string>  GFF3 file annotation file                ($gff_in)
 -fa       <string>  (multi) fasta file holding sequence      ($fa_in)
 -offset   <int>     offset for determination of frame        ($offset)
 -inframe            print only motifs in current ORF         ($inframe)
EOF
exit;
}
