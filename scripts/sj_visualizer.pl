#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-01 23:30:45 mtw>
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
use IPC::Cmd qw(can_run);
use Bio::ViennaNGS  qw(bed2bigBed);
use Bio::ViennaNGS::SpliceJunc qw(bed6_ss_to_bed12);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my ($path_out,$path_annot,$path_ss,$fastaO);
my $mincov = 10;
my $window = 0;
my $outdir = "./";
my $want_bigbed = 0;
my $want_circular=0;
my $s_in = '-';
my $cs_in = '-';
my @result = ();

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("s=s"     => \$s_in,
					   "c=s"     => \$cs_in,
					   "b"       => sub{$want_bigbed=1},
					   "o=s"     => \$outdir,
					   "r=s"     => \$mincov,
					   "w=s"     => \$window,
					   "z"       => sub{$want_circular=1},
					   "man"     => sub{pod2usage(-verbose => 2)},
					   "help|h"  => sub{pod2usage(1)}
					  );


unless ($s_in =~ /^\//) {$s_in = "./".$s_in;}
unless (-f $s_in){
  warn "Could not find input file $s_in given via -s option";
  pod2usage(-verbose => 0);
}

if($want_bigbed==1){
  unless ($cs_in =~ /^\//) {$cs_in = "./".$cs_in;}
  unless (-f $cs_in){
    warn "Could not find input file $cs_in given via -c option";
    pod2usage(-verbose => 0);
  }
  my $bedTobigBed = can_run('bedToBigBed') or
    die 'bedToBigBed utility not found!';
}

#TODO check if we are allowed to write to $outdir
unless ($outdir =~ /\/$/){$outdir .= "/";}
unless (-d $outdir){mkdir $outdir or die $!;}

# Make BED12 line from each splice junction
@result = bed6_ss_to_bed12($s_in,$outdir,$window,$mincov,$want_circular);

my ($bed12) = @result;

if($want_bigbed){
  my $bb12 = bed2bigBed($bed12,$cs_in,$outdir,undef);
}

__END__


=head1 NAME

sj_visualizer.pl - Produce BED12 from BED6 splice junction files.

=head1 SYNOPSIS

sj_visualizer.pl [-s I<FILE>] [-f I<FILE>]
[options]

=head1 DESCRIPTION

Convert splice junctions from mapped RNA-seq data to BED12 format for
easy visualization in genome browsers. The program expects BED6 input
in segemehl splice junction format (see the L<segemehl
documentation|http://www.bioinf.uni-leipzig.de/Software/segemehl/segemehl_manual_0_1_7.pdf>
for details).

=head1 OPTIONS

=over

=item B<-s>

Splice junctions from mapped RNA-seq data in BED6 format

=item B<-c>

Chromosome sizes files

=item B<-b>

Convert resuting BED6 files to bigBed format

=item B<-r>

Minimum read coverage for a splice junction. Only splice junctions
that are supported by at least this number of reads are considered.

=item B<-w>

Expand splice junctions by a window of this size on both sides

=item B<-z>

Filter circular splice junctions

=item B<-o>

Relative output path


=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut

