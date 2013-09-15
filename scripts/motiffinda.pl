#!/usr/bin/perl -w
# -*-CPerl-*-
# Last changed Time-stamp: <2013-09-16 00:27:38 mtw>
#
# Find motifs in annotated sequences. The motif is provided as regex
# via the command line
#
# usage: motiffinda.pl -motif <REXGEX> -gff <GFFFILE> -fa <FASTAFILE>
#                      -inframe
#
# ***********************************************************************
# *  Copyright notice
# *
# *  Copyright 2013 Michael Thomas Wolfinger <michael@wolfinger.eu>
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
use POSIX qw(strftime);
use local::lib;
use ViennaNGS;
use ViennaNGS::AnnoC;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my ($gff_in,$fa_in);
my $motif      = undef;
my $inframe    = 0;
my $timestamp  = strftime("%Y%m%d-%H%M", localtime(time));  # timestamp
my @fa_ids     = ();

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

Getopt::Long::config('no_ignore_case');
&usage() unless GetOptions("gff=s"      => \$gff_in,
                           "fa=s"       => \$fa_in,
                           "motif=s"    => \$motif,
			   "inframe"    => \$inframe,
                           "-help"      => \&usage,
                           "v");
unless ($gff_in =~ /^\//) {$gff_in = "./".$gff_in;}
die "$gff_in not found\n" unless (-f $gff_in);

unless ($fa_in =~ /^\//) {$fa_in = "./".$fa_in;}
die "$fa_in not found\n" unless (-f $fa_in);

unless (defined $motif) {die "please provide motif as regex\n";}

# 1) parse GFF annotation
parse_gff($gff_in);
#print Dumper($feat);
#print Dumper($fstat);

# 2) parse Fasta file, populate @fastaids
get_fasta_dbobjects($fa_in);
#print Dumper(\@fastaids);
foreach my $id (@fastaids) {
  my $obj = $fastadb->get_Seq_by_id($id);  # Bio::PrimarySeq::Fasta
  print Dumper($obj);
}

# 3) find motif in annotates regions
find_motifs();

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

sub find_motifs {

}


sub usage {
  print <<EOF;
Find  motifs in annotated sequences

usage: $0 [options]
program specific options:                                    default:
 -motif    <string>  motif to search for (as regex)           ($motif)
 -gff      <string>  GFF3 file annotation file                ($gff_in)
 -fa       <string>  (multi) fasta file holding sequence      ($fa_in)
 -inframe            print only motifs in current ORF         ($inframe)
EOF
exit;
}
