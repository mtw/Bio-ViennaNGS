#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2013-12-09 14:05:50 mtw>
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
my ($bam_p,$bam_n,$bed_p,$bed_n);
my ($rev,$wantuniq,$wantbed,$bw) = (0)x4;
my $logfile = "bam_split.log";
my $chromsi = undef;
my $bam     = undef;
my @result = ();

Getopt::Long::config('no_ignore_case');
&usage() unless GetOptions("bam=s"           => \$bam,
			   "bed"             => sub{$wantbed = 1},
			   "bw"              => sub{$bw = 1},
			   "c=s"             => \$chromsi,
			   "r"               => sub{$rev = 1},
			   "u"               => sub{$wantuniq = 1},
			   "log=s"           => \$logfile,
                           "-help"           => \&usage,
                           "v");

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

die "ERROR: no BAM file provided" unless (defined $bam);
if ($bw == 1) {
  die "ERROR: chrom_sizes file needed for generating BigWig coverage profiles\n"
    unless (defined $chromsi);
  unless ($wantbed == 1){
    warn "setting -bed option; BED files required for BigWig conversion will be created\n";
    $wantbed = 1;
  }
}

$logfile = $bam . ".bam_split.log";
@result = split_bam($bam,$rev,$wantuniq,$wantbed,$logfile);
$bam_p = $result[0]; # BAM file containing fragments of [+] strand
$bam_n = $result[1]; # BAM file containing fragments of [-] strand
$bed_p = $result[2]; # BED file containing fragments of [+] strand
$bed_n = $result[3]; # BED file containing fragments of [-] strand

if ($bw == 1) {
  bam2bw($bam_p,$chromsi);
  bam2bw($bam_n,$chromsi);
}

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
 -bed               create BED file for each split BAM      ($wantbed)
 -bw                create BedGraph and BigWig files        ($bw)
 -c                 chrom_sizes for generating BigWigs      ($chromsi)
 -r                 reverse +/- strand mapping (due to      ($rev)
                    RNA-seq configuration)
 -u                 filter unique alignemnts                ($wantuniq)
 -log      <file>   log file                                ($logfile)
 -help               print this information

EOF
exit;
}
