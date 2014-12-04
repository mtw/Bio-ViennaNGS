#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-04 12:52:31 mtw>
#
# Split BAM files according to their strands, optionally filter unique
# mappers
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

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use Bio::ViennaNGS qw(split_bam bed_or_bam2bw);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
my ($bam_p,$bam_n,$bed_p,$bed_n,$size_p,$size_n,$basename,$bamdir,$bamext,$cmd);
my ($rev,$wantuniq,$wantbed,$wantnorm,$bw) = (0)x5;
my $logfile = "bam_split.log";
my $chromsi = undef;
my $fullbam = undef;
my $destdir = "./";
my $scale   = 1000000;
my @result  = ();
my $this_function = (caller(0))[3];

Getopt::Long::config('no_ignore_case');
&usage() unless GetOptions("bam=s"           => \$fullbam,
			   "bed"             => sub{$wantbed = 1},
			   "bw"              => sub{$bw = 1},
			   "c=s"             => \$chromsi,
			   "norm"            => sub{$wantnorm = 1},
			   "o=s"             => \$destdir,
			   "r"               => sub{$rev = 1},
			   "s=s"             => \$scale,
			   "u"               => sub{$wantuniq = 1},
			   "log=s"           => \$logfile,
                           "-help"           => \&usage,
                           "v");

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

die "ERROR [$this_function] No BAM file provided" unless (defined $fullbam);
if ($bw == 1) {
  die "ERROR [$this_function] chrom_sizes file needed for generating BigWig coverage profiles\n"
    unless (defined $chromsi);
  unless ($wantbed == 1){$wantbed = 1;}
}
unless ($destdir =~ /\/$/){$destdir .= "/";}
unless (-d $destdir){$cmd = "mkdir -p $destdir"; system($cmd);}
unless ($fullbam =~ /^\// || $fullbam =~ /\.\//){$fullbam = "./".$fullbam;}
($basename,$bamdir,$bamext) = fileparse($fullbam,qr/\.[^.]*/);

$logfile = $destdir.$basename.".bam_split.log";
@result = split_bam($fullbam,$rev,$wantuniq,$wantbed,$destdir,$logfile);
$bam_p  = $result[0]; # BAM file containing fragments of [+] strand
$bam_n  = $result[1]; # BAM file containing fragments of [-] strand
$size_p = $result[2]; # of alignments on [+] strand
$size_n = $result[3]; # of alignments on [-] srand
$bed_p  = $result[4]; # BED file containing fragments of [+] strand
$bed_n  = $result[5]; # BED file containing fragments of [-] strand

if ($bw == 1) {
  $destdir = $destdir."vis";
  $cmd = "mkdir -p $destdir"; system($cmd);
  bed_or_bam2bw("bed",$bed_p,$chromsi,"+",$destdir,$wantnorm,$size_p,$scale,$logfile);
  bed_or_bam2bw("bed",$bed_n,$chromsi,"-",$destdir,$wantnorm,$size_n,$scale,$logfile);
}

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

sub usage {
 print <<EOF;

bam_split.pl:  Split a BAM file according to strands.

Optionally filter unique alignments by inspecting NH:i SAM attribute
Optionally create bedGraph and (stranded |normalized) bigWig coverage
for UCSC visualization

usage: $0 -bam <BAMFILE> [options]
program specific options:                                   default:
 -bam   <file>   specify BAM file                           ($fullbam)
 -bed            create BED file for each split BAM         ($wantbed)
 -bw             create BedGraph and bigWig files           ($bw)
 -c     <file>   chrom_sizes for generating bigWig files    ($chromsi)
 -norm           normalize resulting bigWig files           ($wantnorm)
 -o     <path>   output directory                           ($destdir)
 -r              reverse +/- strand mapping (according      ($rev)
                 to RNA-seq library preparation protocol)
 -s     <int>    scale bigWig files to this number          ($scale)
 -u              filter unique alignemnts                   ($wantuniq)
 -log   <file>   log file                                   ($logfile)
 -help           print this information

EOF
exit;
}
