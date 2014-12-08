#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-09 00:45:17 mtw>
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
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Pod::Usage;
use Data::Dumper;
use File::Basename;
use Path::Class;
use Bio::ViennaNGS qw(split_bam bed_or_bam2bw);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
my ($bam_p,$bam_n,$bed_p,$bed_n,$size_p,$size_n,$basename,$bamdir,$bamext,$lf);
my ($rev,$wantuniq,$wantbed,$wantnorm,$bw) = (0)x5;
my $logext = ".bam_split.log";
my $cs_in = "-";
my $bam_in = "-";
my $outdir = "./";
my $scale   = 1000000;
my @result  = ();

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("bam=s"      => \$bam_in,
					   "bed"        => sub{$wantbed = 1},
					   "bw"         => sub{$bw = 1},
					   "cs=s"       => \$cs_in,
					   "norm"       => sub{$wantnorm = 1},
					   "o|out=s"    => \$outdir,
					   "r|reverse"  => sub{$rev = 1},
					   "scale=s"    => \$scale,
					   "uniq"       => sub{$wantuniq = 1},
					   "l|log=s"    => \$logext,
					   "man"        => sub{pod2usage(-verbose => 2)},
					   "help|h"     => sub{pod2usage(1)}
					  );

unless ($bam_in =~ /^\// || $bam_in =~ /\.\//){$bam_in = "./".$bam_in;}
unless (-f $bam_in){
  warn "Could not find input file $bam_in given via --bam option";
  pod2usage(-verbose => 0);
}

if ($bw == 1) {
  unless (-f $cs_in){
    warn "Could not find input file $cs_in given via --cs option";
    pod2usage(-verbose => 0);
  }
  unless ($wantbed == 1){$wantbed = 1;}
}

#TODO check if we are allowed to write to $outdir
#unless ($outdir =~ /\/$/){$outdir .= "/";}
unless (-d $outdir){mkdir $outdir or die $!;}

($basename,$bamdir,$bamext) = fileparse($bam_in,qr/\..*/);

$lf = file($outdir,$basename.$logext);

@result = split_bam($bam_in,$rev,$wantuniq,$wantbed,$outdir,$lf);
$bam_p  = $result[0]; # BAM file containing fragments of [+] strand
$bam_n  = $result[1]; # BAM file containing fragments of [-] strand
$size_p = $result[2]; # of alignments on [+] strand
$size_n = $result[3]; # of alignments on [-] srand
$bed_p  = $result[4]; # BED file containing fragments of [+] strand
$bed_n  = $result[5]; # BED file containing fragments of [-] strand

if ($bw == 1) {
  my $od = dir($outdir,"vis");
  mkdir $od or die $!;
  bed_or_bam2bw("bed",$bed_p,$cs_in,"+",$od,$wantnorm,$size_p,$scale,$lf);
  bed_or_bam2bw("bed",$bed_n,$cs_in,"-",$od,$wantnorm,$size_n,$scale,$lf);
}


__END__


=head1 NAME

bam_split.pl - Split a BAM file by strands

=head1 SYNOPSIS

bam_split.pl [--bam I<FILE>] [options]

=head1 DESCRIPTION

Split a BAM file by strands and create two new BAM file: One
containing all reads that map to the positive strand and another one
with all reads mapped to the negative strand. Optionally filter unique
alignments by inspecting NH:i SAM attribute.

Optionally create bedGraph and (stranded |normalized) bigWig coverage
for UCSC visualization

=head1 OPTIONS

=over

=item B<--bam>

Input file in BAM format

=item B<--bed>

Create a BED6 file for each split BAM file

=item B<--bw>

Create BedGraph and bigWig coverage files for e.g. genome browser
visualization.

=item B<--cs>

Chromosome sizes file (required if B<--bw> is given).

=item B<--norm>

Normalize resulting bigWig files

=item B<--out -o>

Output directory

=item B<--reverse -r>

Reverse the +/- strand mapping. This is required to achieve proper
strand assignments for certain RNA-seq library preparation protocol.

=item B<--scale>

If B<--bw> is given, scale bigWig files to this number. Default is 1000000.

=item B<--uniq>

Filter uniquely mapped reads by inspecting the NH:i: SAM
attribute. See also the I<bam_uniq.pl> utility, which extracts both
uniquely and multiply mapped reads from BAM files without
strand-splitting.

=item B<--log -l>

Log file extension. Default is ".bam_split.log". The log file is
created in the directory given via B<-o> and its name is constructed
from the base name of the input BAM file and the log filename extension.

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut

