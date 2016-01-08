#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2016-01-08 14:04:59 mtw>
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
use Path::Class;
use File::Basename;
use Bio::ViennaNGS::Bam qw(uniquify_bam);
use Bio::ViennaNGS::Util  qw(mkdircheck);
use Cwd;
use Data::Dumper;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my ($bam_in,$cwd,$lf,$basename,$bamdir,$bamext);
my $logext = ".bam_uniq.log";
my $outdir = "./";
my @result = ();

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("bam=s"     => \$bam_in,
					   "o|out=s"   => \$outdir,
					   "l|log=s"   => \$logext,
					   "man"       => sub{pod2usage(-verbose => 2)},
					   "help|h"    => sub{pod2usage(1)}
					  );

unless (-f $bam_in){
  warn "Could not find input file $bam_in given via -bam option";
  pod2usage(-verbose => 0);
}
$cwd = getcwd();
unless ($bam_in =~ /^\// || $bam_in =~ /\.\//) {$bam_in = file($cwd,$bam_in);}
unless (-d $outdir){mkdircheck($outdir)};

($basename,$bamdir,$bamext) = fileparse($bam_in,qr/\.bam/);
$lf = file($outdir,$basename.$logext);

@result = uniquify_bam($bam_in,$outdir,$lf);

__END__


=head1 NAME

bam_uniq.pl - Deconvolute BAM files into unique and multi mappers

=head1 SYNOPSIS

bam_uniq.pl [-bam I<FILE>] [-o I<DIR>]
[options]

=head1 DESCRIPTION

This program extracts uniquely and multiply mapping reads from BAM
files acoording to the NH:i SAM attribute and creates separate BAM
files for unique and multi mappers, respectively.

=head1 OPTIONS

=over

=item B<--bam>

BAM file to extract unique and multi mappers from

=item B<--out -o>

Output path

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut

