#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2015-10-23 12:36:14 mtw>
#
# Find peaks/enriched regions of certain size in RNA-seq data
#
# usage: rnaseq_peakfinder --bgpos pos.bg --bgneg neg.pg
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
use Bio::ViennaNGS::FeatureIO;
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
					   "l|length=i"      => \$maxlen,
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
my $io_pos = Bio::ViennaNGS::FeatureIO->new(file => "$infile1", # [+] strand
					    filetype => "BedGraph",
					    objecttype => "BedGraph",
					   );

my $io_neg = Bio::ViennaNGS::FeatureIO->new(file => "$infile2", # [-] strand
					    filetype => "BedGraph",
					    objecttype => "BedGraph",
					   );

# get an instance of Bio::ViennaNGS::Peak
my $peaks = Bio::ViennaNGS::Peak->new(winsize    => $winsize,
				      interval   => $wininterval,
				      mincov     => $mincov,
				      length     => $maxlen,
				      threshold  => $threshold,
				     );

# parse BedGraph coverage data from FeatureIO into the Peak object
$peaks->populate_data($io_pos, $io_neg);

# identify covered regions
$peaks->raw_peaks($dest, "rnaseq_peakfinder", $lf);

# characterize peaks
$peaks->final_peaks($dest, "rnaseq_peakfinder", $lf);


__END__


=head1 NAME

rnaseq_peakfinder.pl - Identify peaks/enriched regions in RNA-seq data

=head1 SYNOPSIS

rnaseq_peakfinder.pl [--bgpos I<FILE>] [--bgneg I<FILE>] [options]

=head1 DESCRIPTION

This program identifies peaks in RNA-seq data. Starting from coverage
information in bedGraph format, this tools applies a two-step sliding
window approach to characterize enriched regions with predefined
properties, including maximum length, minimum coverage or maximum
coverage at both enads of the genomic inerval.

B<Please note>: It is highly recommended to use I<normalized> input
data.

=head1 OPTIONS

=over

=item B<--bgpos>

BedGraph input file containing coverage of the [+] strand.

=item B<--bgneg>

BedGraph input file containing coverage of the [-] strand.

=item B<--winsize -w>

Size of the sliding window in nt.

=item B<--interval -i>

Size of the interval the sliding window is shifted at each step ('step
size').

=item B<--mincov -m>

Minimum coverage required for an enriched region to be considered.

=item B<--length -l>

Maximum length of a peak in nt.

=item B<--threshold -t>

Percentage of the maximum coverage value allowed at both ends of the
peaks (default 0.1). This value is used to identify peak boundaries.

=item B<--out -o>

Output directory.

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 NOTES

The memory footprint of this tool is rather high (several GB for
eucaryotic systems). This is due to the fact that the input BedGraph
files are first parsed into an Array of
L<Bio::ViennaNGS::BedGraphEntry> objects by
L<Bio::ViennaNGS::FeatureIO>. In a second step, this array is parsed
into a Hash of Arrays data structure within L<Bio::ViennaNGS::Peaks>
to allow for efficient window sliding. This may be refactored in a
future release.


=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut





