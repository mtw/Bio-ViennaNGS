#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-09-24 16:12:50 mtw>
#
# ***********************************************************************
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
use local::lib;
#use ViennaNGS;
#use ViennaNGS::AnnoC;
use ViennaNGS::SpliceJunc;
#use Bio::Seq;
#use Bio::SeqIO;


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my ($path_out,$path_annot,$path_ss);
my $mincov = 10;
my $window = 10;
my $max_intron_length = 7000; # 7000 for A. thaliana
my $destdir = "./";
my $prefix = "";
my $dirname_annot = "annotated";
my $dirname_ss = "mapped";
my $want_canonical = 0;
my $s_in = '-';
my $bed12_in = '-';

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("s=s"     => \$s_in,
					   "a=s"     => \$bed12_in,
					   "c"       => sub{$want_canonical=1},
					   "d=s"     => \$destdir,
					   "r=s"     => \$mincov,
					   "i=s"     => \$max_intron_length,
					   "p=s"     => \$prefix,
					   "w=s"     => \$window,
					   "man"     => sub{pod2usage(-verbose => 2)},
					   "help|h"  => sub{pod2usage(1)}
					  );

unless ($bed12_in =~ /^\//) {$bed12_in = "./".$bed12_in;}
unless (-f $bed12_in) {
  warn "Could not find input file $bed12_in given via -a option";
  pod2usage(-verbose => 0);
}

unless ($s_in =~ /^\//) {$s_in = "./".$s_in;}
unless (-f $s_in){
  warn "Could not find input file $s_in given via -s option";
  pod2usage(-verbose => 0);
}

unless ($destdir =~ /\/$/){$destdir .= "/";}
unless (-d $destdir){my $cmd = "mkdir -p $destdir"; system($cmd);}
$path_annot = $destdir.$dirname_annot;
unless (-d $path_annot){my $cmd = "mkdir -p $path_annot"; system($cmd);}
$path_ss = $destdir.$dirname_ss;
unless (-d $path_ss){my $cmd = "mkdir -p $path_ss"; system($cmd);}

# Extract annotated splice sites from a bed12 file:
# Process annotation and create bed6 files of annotated splice
# sites (+/- 10nt) FOR EACH TRANSCRIPT. Name the files like
# chr1-3630-5899-AT1G010101.1.bed.
# To get those files, process either the genePred file or explicitly
# create the intronic regions via 'bed12ToBed6 -i tr1.bed >
# tr1_exons.bed' and 'subtractBed -a tr1.bed -b tr1_exons.bed -s'
bed6_ss_from_bed12($bed12_in,$path_annot,$window);

# Extract mapped splice junctions and create a BED6 file for each
# of them
bed6_ss_from_rnaseq($s_in,$path_ss,$window,$mincov);

# Check for each splice junction seen in RNA-seq if it overlaps with
# any annotated splice junction
intersect_sj($path_annot,$path_ss,$destdir,$prefix,$window,$max_intron_length,$want_canonical);


__END__


=head1 NAME

splice_site_summary.pl - Find novel splice junctions in RNA-seq data.

=head1 SYNOPSIS

splice_site_summary.pl [-s I<FILE>] [-a I<FILE>] [options] [-man] [-help]

=head1 OPTIONS

=over 

=item B<-a>

Reference annotation as BED12

=item B<-s>

Splice junctions from mapped RNA-seq data as BED6

=item B<-c>

Filter I<canonical> splice junctuions

=item B<-i>

Maximum intron length. Splice junctions covering introns larger than
this value are not considered.

=item B<-r>

Minimum read coverage for a splice junction. Only splice junctions
that are supported by at least this number of reads are considered.

=item B<-w>

Expand splice junctions by a window of this size on both sides

=item B<-d>

Output path

=item B<-p>

Prefix all output files with this string

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 DESCRIPTION

This program identifies and characterizes splice sites from mapped
RNA-seq data against annotated splice junctions.

=head1 AUTHORS

Michael Thomas Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut

