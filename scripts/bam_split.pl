#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-04-04 10:10:39 mtw>
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
use File::Basename;
use ViennaNGS;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
my ($bam_p,$bam_n,$bed_p,$bed_n,$basename,$bamdir,$bamext,$cmd);
my ($rev,$wantuniq,$wantbed,$bw) = (0)x4;
my $logfile = "bam_split.log";
my $chromsi = undef;
my $fullbam = undef;
my $destdir = "./";
my @result = ();

Getopt::Long::config('no_ignore_case');
&usage() unless GetOptions("bam=s"           => \$fullbam,
			   "bed"             => sub{$wantbed = 1},
			   "bw"              => sub{$bw = 1},
			   "c=s"             => \$chromsi,
			   "o=s"             => \$destdir,
			   "r"               => sub{$rev = 1},
			   "u"               => sub{$wantuniq = 1},
			   "log=s"           => \$logfile,
                           "-help"           => \&usage,
                           "v");

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

die "ERROR: no BAM file provided" unless (defined $fullbam);
if ($bw == 1) {
  die "ERROR: chrom_sizes file needed for generating BigWig coverage profiles\n"
    unless (defined $chromsi);
  unless ($wantbed == 1){
    warn "setting -bed option; BED files required for BigWig conversion will be created\n";
    $wantbed = 1;
  }
}
unless ($destdir =~ /\/$/){$destdir .= "/";}
unless (-d $destdir){$cmd = "mkdir -p $destdir"; system($cmd);}
unless ($fullbam =~ /^\// || $fullbam =~ /\.\//){$fullbam = "./".$fullbam;}
($basename,$bamdir,$bamext) = fileparse($fullbam,qr/\.[^.]*/);

$logfile = $destdir.$basename.".bam_split.log";
@result = split_bam($fullbam,$rev,$wantuniq,$wantbed,$destdir,$logfile);
$bam_p = $result[0]; # BAM file containing fragments of [+] strand
$bam_n = $result[1]; # BAM file containing fragments of [-] strand
$bed_p = $result[2]; # BED file containing fragments of [+] strand
$bed_n = $result[3]; # BED file containing fragments of [-] strand

if ($bw == 1) {
  $destdir = $destdir."vis";
  $cmd = "mkdir -p $destdir"; system($cmd);
  bed2bw($bed_p,$chromsi,"+",$destdir);
  bed2bw($bed_n,$chromsi,"-",$destdir);
}

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

sub usage {
 print <<EOF;

bam_split.pl:  Split a BAM file according to strands.

Optionally filter unique alignments by inspecting NH:i SAM attribute
Optionally create BedGraph and stranded BigWig coverage for UCSC visualization

usage: $0 -bam <BAMFILE> [options]
program specific options:                                   default:
 -bam   <file>   specify BAM file                           ($fullbam)
 -bed            create BED file for each split BAM         ($wantbed)
 -bw             create BedGraph and BigWig files           ($bw)
 -c     <file>   chrom_sizes for generating BigWigs         ($chromsi)
 -o     <path>   output directory                           ($destdir)
 -r              reverse +/- strand mapping (according      ($rev)
                 to RNA-seq library preparation protocol)
 -u              filter unique alignemnts                   ($wantuniq)
 -log   <file>   log file                                   ($logfile)
 -help           print this information

EOF
exit;
}
