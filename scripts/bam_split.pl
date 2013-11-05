#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2013-11-05 15:02:42 mtw>
#
# Split BAM files according to their strands, optionally filter unique mappers
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

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use ViennaNGS;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
my ($bam_pos,$bam_neg);
my ($rev,$uniq,$bw) = (0)x3;
my $logfile = "bam_split.log";
my $chromsi = undef;
my $bam     = undef;

Getopt::Long::config('no_ignore_case');
&usage() unless GetOptions("bam=s"           => \$bam,
			   "bw"              => sub{$bw = 1},
			   "c=s"             => \$chromsi,
			   "r"               => sub{$rev = 1},
			   "u"               => sub{$uniq = 1},
			   "log=s"           => \$logfile,
                           "-help"           => \&usage,
                           "v");

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

die "ERROR: no BAM file provided" unless (defined $bam);
die "ERROR: chrom_sizes file needed for generating BigWig coverage profiles\n"
  unless (defined $chromsi && $bw == 1);

$logfile = $bam . ".bam_split.log";
($bam_pos,$bam_neg) = split_bam($bam,$rev,$uniq,$logfile);
bam2bw($bam_pos,$chromsi);
bam2bw($bam_neg,$chromsi);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

sub usage {
 print <<EOF;

bam_split.pl:  Split a BAM file according to strands.

Optionally filter unique alignments by inspecting NH:i SAM attribute
Optionally create BedGraph and BigWig coverage for UCSC visualization

usage: $0 -bam <BAMFILE> [options]
program specific options:                                   default:
 -bam      <string> specify BAM file                        ($bam)
 -bw                create BedGraph and BigWig files        ($bw)
 -c                 chrom_sizes for generating BigWigs      ($chromsi)
 -r                 reverse +/- strand mapping (due to      ($rev)
                    RNA-seq configuration)
 -u                 filter unique alignemnts                ($uniq)
 -log      <file>   log file                                ($logfile)
 -help               print this information

EOF
exit;
}
