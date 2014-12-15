#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-13 01:05:02 mtw>
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
use Bio::ViennaNGS::Util  qw(bed2bigBed);
use Bio::ViennaNGS::SpliceJunc qw(bed6_ss_from_bed12 bed6_ss_from_rnaseq intersect_sj);
use Bio::ViennaNGS::Fasta;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my ($path_out,$path_annot,$path_ss,$fastaO);
my $mincov = 10;
my $window = 10;
my $max_intron_length = 7000; # 7000 for A. thaliana
my $outdir = "./";
my $prefix = "";
my $dirname_annot = "annotated";
my $dirname_ss = "mapped";
my $want_canonical = 0;
my $want_bigbed = 0;
my $s_in = '-';
my $bed12_in = '-';
my $fa_in = '-';
my $cs_in = '-';
my %fastaobj = ();
my @result = ();

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("s=s"     => \$s_in,
					   "a=s"     => \$bed12_in,
					   "f=s"     => \$fa_in,
					   "c=s"     => \$cs_in,
					   "n"       => sub{$want_canonical=1},
					   "b"       => sub{$want_bigbed=1},
					   "o=s"     => \$outdir,
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

if ($want_canonical==1){
  unless ($fa_in =~ /^\//) {$fa_in = "./".$fa_in;}
  unless (-f $fa_in){
    warn "Could not find input file $fa_in given via -f option";
    pod2usage(-verbose => 0);
  }
}

if($want_bigbed==1){
  unless ($cs_in =~ /^\//) {$cs_in = "./".$cs_in;}
  unless (-f $cs_in){
    warn "Could not find input file $cs_in given via -c option";
    pod2usage(-verbose => 0);
  }
}

#TODO check if we are allowed to write to $outdir
unless ($outdir =~ /\/$/){$outdir .= "/";}
unless (-d $outdir){mkdir $outdir or die $!;}
$path_annot = $outdir.$dirname_annot;
unless (-d $path_annot){mkdir $path_annot or die $!;}
$path_ss = $outdir.$dirname_ss;
unless (-d $path_ss){mkdir $path_ss or die $!;}

if($want_canonical){
  $fastaO = Bio::ViennaNGS::Fasta->new(fa=>$fa_in);
}

# Extract annotated splice sites from a BED12 file
bed6_ss_from_bed12($bed12_in,$path_annot,$window,$want_canonical,$fastaO);

# Extract mapped splice junctions and create a BED6 file for each
# of them
bed6_ss_from_rnaseq($s_in,$path_ss,$window,$mincov,$want_canonical,$fastaO);

# Check for each splice junction seen in RNA-seq if it overlaps with
# any annotated splice junction
@result = intersect_sj($path_annot,$path_ss,$outdir,$prefix,$window,$max_intron_length);
my ($exist,$novel) = @result;

if($want_bigbed){
  my $bbe = bed2bigBed($exist,$cs_in,$outdir,undef);
  my $bbn = bed2bigBed($novel,$cs_in,$outdir,undef);
}

__END__


=head1 NAME

splice_site_summary.pl - Find novel splice junctions in RNA-seq data.

=head1 SYNOPSIS

splice_site_summary.pl [-s I<FILE>] [-a I<FILE>] [-f I<FILE>]
[options]

=head1 DESCRIPTION

This program identifies and characterizes splice sites from mapped
RNA-seq data against annotated splice junctions.

=head1 OPTIONS

=over

=item B<-a>

Reference annotation in BED12 format

=item B<-s>

Splice junctions from mapped RNA-seq data in BED6 format

=item B<-f>

Reference genome in Fasta format

=item B<-c>

Chromosome sizes files

=item B<-b>

Convert resuting BED6 files to bigBed format

=item B<-n>

Filter I<canonical> splice junctions

=item B<-i>

Maximum intron length. Splice junctions covering introns larger than
this value are not considered.

=item B<-r>

Minimum read coverage for a splice junction. Only splice junctions
that are supported by at least this number of reads are considered.

=item B<-w>

Expand splice junctions by a window of this size on both sides

=item B<-o>

Relative output path

=item B<-p>

Prefix all output files with this string

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut

