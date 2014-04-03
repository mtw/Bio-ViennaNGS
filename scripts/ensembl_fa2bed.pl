#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2013-11-05 16:06:38 mtw>
#
# Extract annotation from ENSEMBL FASTA files and convert to BED
#
# usage: ensembl_fa2bed.pl -fa <FASTAFILE>
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

my ($gff_in,$fa_in,$obj);
my $motif      = undef;
my $inframe    = 0;
my $offset     = 0;
my $timestamp  = strftime("%Y%m%d-%H%M", localtime(time));  # timestamp
my @fa_ids     = ();

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

Getopt::Long::config('no_ignore_case');
&usage() unless GetOptions(
                           "fa=s"       => \$fa_in,
                           "-help"      => \&usage,
                           "v");

unless ($fa_in =~ /^\//) {$fa_in = "./".$fa_in;}
die "$fa_in not found\n" unless (-f $fa_in);

# 1) parse Fasta file, populate @fastaids
get_fasta_ids($fa_in);
#print Dumper(\@fastaids);
foreach my $id (@fastaids) {

  my $hd = $fastadb->header($id);
  print Dumper($hd);
 $obj = $fastadb->get_Seq_by_id($id);  # Bio::PrimarySeq::Fasta
 print Dumper($obj);
}

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

sub usage {
  print <<EOF;
Extract annotation information from ENSEMBL FASTA and convert to BED

usage: $0 [options]
program specific options:                                    default:
 -fa       <string>  (multi) fasta file holding sequence      ($fa_in)
EOF
exit;
}
