#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-10 00:34:41 mtw>
#
# Convert GFF3 to BED12; produce separate BED files for each gbkey
# (CDS/tRNA/rRNA/ncRNA etc.)
#
# usage: gff2bed.pl --gff input.gff
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
#

use strict;
use warnings;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Pod::Usage;
use Data::Dumper;
use File::Basename;
use Path::Class;
use Bio::ViennaNGS::AnnoC;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my $feature       = undef; # look for all features per default
my $gff_in        = '-';
my $outdir       = './';
my ($basename,$gffdir,$gffext,$obj);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("gff=s"       => \$gff_in,
					   "o|out=s"     => \$outdir,
					   "f|feature=s" => \$feature,
					   "man"         => sub{pod2usage(-verbose => 2)},
					   "help|h"      => sub{pod2usage(1)}
					  );

unless ($gff_in =~ /^\// || $gff_in =~ /\.\//){$gff_in = "./".$gff_in;}
unless (-f $gff_in){
  warn "Could not find input file $gff_in given via --gff option";
  pod2usage(-verbose => 0);
}

#TODO check if we are allowed to write to $outdir
unless ($outdir =~ /\/$/){$outdir .= "/";}
unless (-d $outdir){mkdir $outdir or die $!;}

($basename,$gffdir,$gffext) = fileparse($gff_in,qr/\..*/);

$obj = Bio::ViennaNGS::AnnoC->new();
$obj->parse_gff($gff_in);
$obj->featstat;
$obj->feature_summary($outdir);
$obj->features2bed($feature,$outdir,$basename,undef);


__END__


=head1 NAME

gff2bed.pl - Convert (non-spliced) GFF3 to BED12

=head1 SYNOPSIS

gff2bed.pl [--gff I<FILE>] [options]

=head1 DESCRIPTION

Convert feature annotation of non-spliced organisms in GFF3 format to
BED12. A separate BED12 file will be created for each genomic feature
type (eg CDS, tRNA, rRNA, ncRNA, etc.).

This script serves as a reference implementation of code fragments
from bio::ViennaNGS::AnnoC and has been successfully tested with
NCBI's bacterial annotation in GFF3 format.

Please note that this script is NOT meant to be run on GFF genome
annotation other than bacteria since it DOES NOT consider splice
events present in most higher organisms.

=head1 OPTIONS

=over

=item B<--gff>

Input GFF file.

=item B<--feature -f>

Specify feature type (eg. CDS,tRNA,rRNA,SBS, etc) to be extracted from GFF3.

=item B<--out -o>

Output path.

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut


