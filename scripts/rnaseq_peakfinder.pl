#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2015-10-14 16:10:02 mtw>
#
# Find peaks/enriched regions of certain size in RNA-seq data
#
# ***********************************************************************
# *  Copyright notice
# *
# *  Copyright 2015 Michael T. Wolfinger <michael@wolfinger.eu>
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
use Cwd;
use Path::Class;
use Bio::ViennaNGS::Peak;


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
my ($infile1, $infile2, $lf);
my $logname            = "rnaseq_peakfinder.log.txt";
my $fn_rawpeaks        = "rawpeaks.bed";
my $fn_candidatepeaks  = "candidatepeaks.bed";
my $winsize            = 15;           # sliding window size
my $wininterval        = 1;            # windows sliding interval
my $mincov             = 100;          # minimum coverage
my $maxlen             = 100;          # maximum peak length
my $threshold          = 0.1;          # cut off threhsold for peak end detection
my $debug              = 0;            # debug switch
my $dest               = getcwd();


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("bgpos=s"         => \$infile1,
					   "bgneg=s"         => \$infile2,
					   "w|winsize=i"     => \$winsize,
					   "i|interval=i"    => \$wininterval,
					   "m|mincov=i"      => \$mincov,
					   "e|length=i"      => \$maxlen,
					   "t|threshold=f"   => \$threshold,
					   "out|o=s"         => \$dest,
					   "d|debug"         => sub{$debug=1},
                                           "man"             => sub{pod2usage(-verbose => 2)},
                                           "help|h"          => sub{pod2usage(1)}
                                          );

unless(-f $infile1){
  warn "Could not find input bedGraph file for [+] strand provided via ---bgpos option";
  pod2usage(-verbose => 0);
}
unless(-f $infile2){
  warn "Could not find input bedGraph file for [-] strand provided via ---bgneg option";
  pod2usage(-verbose => 0);
}
my $cwd = getcwd();
unless ($infile1 =~ /^\// || $infile1 =~ /^\.\//){$infile1 = file($cwd,$infile1);}
unless ($infile2 =~ /^\// || $infile2 =~ /^\.\//){$infile2 = file($cwd,$infile2);}
unless ($dest =~ /\/$/){$dest .= "/";}
unless (-d $dest){mkdir $dest;}

$lf = file($dest,$logname);

# parse input files
my $peaks = Bio::ViennaNGS::Peak->new(winsize    => $winsize,
				      interval   => $wininterval,
				      mincov     => $mincov,
				      length     => $maxlen,
				      threshold  => $threshold,
				     );
$peaks->parse_coverage_bedgraph($infile1,$infile2);

$peaks->raw_peaks($dest, "rnaseq_peakfinder", $lf);

#print Dumper(\$peaks);

