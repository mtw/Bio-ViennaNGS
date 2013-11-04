#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2013-11-04 21:42:36 mtw>
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
use local::lib;
use Getopt::Long;
use Data::Dumper;
use ViennaNGS;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
my $bam;
my ($rev,$uniq) = (0)x2;
my $logfile = "bam_split.log";

Getopt::Long::config('no_ignore_case');
&usage() unless GetOptions("bam=s"           => \$bam,
			   "r"               => sub{$rev = 1},
			   "u"               => sub{$uniq = 1},
			   "log=s"           => \$logfile,
                           "-help"           => \&usage,
                           "v");

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

split_bam($bam,$rev,$uniq,$logfile);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#



sub usage {
  print <<EOF;
Split a BAM file according to strands, optionally report unique alignments

usage: $0 -bam <BAMFILE> [options]
program specific options:                                   default:
 -bam          <string> specify BAM file                    ($bam)
 -r                     reverse +/- strand mapping (due to  ($rev)
                        RNA-seq configuration)
 -u                     filter unique alignemnts            ($uniq)
 -log          <file>   log file                            ($logfile)
 -help                  print this information

EOF
exit;
}
