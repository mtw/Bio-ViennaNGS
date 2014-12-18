#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-13 01:01:55 mtw>
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
use Path::Class;
use Bio::ViennaNGS::Util qw(bed_or_bam2bw);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my ($bam_in,$lf,$bwfile);
my $logfile = "bam_to_bigwig.log";
my $outdir = "./";
my $want_bigbed = 0;
my $cs_in = "-";
my $strand = "+";

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("b|bam=s"    => \$bam_in,
					   "c|cs=s"     => \$cs_in,
					   "o=s"        => \$outdir,
					   "l|log=s"    => \$logfile,
					   "s|strand=s" => \$strand,
					   "man"        => sub{pod2usage(-verbose => 2)},
					   "help|h"     => sub{pod2usage(1)}
					  );


unless ($bam_in =~ /^\//) {$bam_in = "./".$bam_in;}
unless (-f $bam_in){
  warn "Could not find input file $bam_in given via --bam option";
  pod2usage(-verbose => 0);
}
unless (-f $cs_in){
  warn "Could not find input file $cs_in given via --cs option";
  pod2usage(-verbose => 0);
}
unless ($strand =~ /[\+\-]/) { 
  warn "Invalid value '$strand' given via -s option. Please specify
  either '+' or '-' for positive or negative strand, respectively";
  pod2usage(-verbose => 0);
 }

#TODO check if we are allowed to write to $outdir
unless ($outdir =~ /\/$/){$outdir .= "/";}
unless (-d $outdir){mkdir $outdir or die $!;}
$lf = file($outdir,$logfile);

$bwfile = bed_or_bam2bw("bam",$bam_in,$cs_in,$strand,$outdir,0,0,1.,$lf);

print "$bwfile\n";

__END__


=head1 NAME

bam_to_bigWig.pl - Make bigWig coverage profiles from BAM files

=head1 SYNOPSIS

bam_to_bigWig.pl [--bam I<FILE>] [--cs I<FILE>] [--strand I<+/->] [options]

=head1 DESCRIPTION

Produce bigWig coverage profiles from (aligned) BAM files, explicitly
considering strandedness. The most natural use case of this tool is to
create strand-aware coverage profiles in bigWig format for genome
browser visualization.

=head1 OPTIONS

=over

=item B<--bam -b>

Input file in BAM format

=item B<--cs -c>

Chromosome sizes file

=item B<--strand -s>

Use this option if the input BAM file is strictly strand-specific,
ie. contains B<only> reads mapped to either the positive or negative
strand. Possible values are either '+' or '-'. If the value given here
is '+', the interim bedGraph file will be created with positive
values. A '-' given here will create the inerim bedGraph file with
negative values, which is required for proper visualization of bigWig
files holding coverage profiles of reads mapped to the negative strand
in the UCSC genome browser. If the input BAM file is not
strand-specific, ie contains reads mapped to both positive and
negative strand, then the default value '+' will be used, resulting in
bigWig coverage profiles rendered in positive (y-axis direction) in
the UCSC genome browser.

=item B<-o>

Output directory

=item B<--log -l>

Name of the log file. Unless specified, the default log file will be
"bam_to_bigwig.log" in the given output directory.


=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut

