#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-10-04 00:16:24 mtw>
#
# Convert GFF3 to BED12; produce separate BED files for each gbkey
# (CDS/tRNA/etc)
#
# usage: gff2bed.pl -i input.gff
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
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use Bio::ViennaNGS::AnnoC qw(&parse_gff &feature_summary &features2bed $fstat $feat);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my $chr_name      = '';
my $feature       = 'gene'; # look for 'gene' per default
my $infile        = '';
my $workdir       = undef;
my $ext           = ".gff";
my $sortBed       = "sortBed";
my ($bn);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

Getopt::Long::config('no_ignore_case');
&usage() unless GetOptions("chr=s"      => \$chr_name,
			   "i=s"        => \$infile,
			   "wd=s"       => \$workdir,
			   "feat=s"     => \$feature,
                           "-help"      => \&usage,
                           "v");
unless (defined $workdir) {
  $workdir = "./";
}
$workdir .= "/" unless ($workdir =~ /\/$/);

$bn = basename($infile, $ext);
$infile = $workdir.$infile;
print "gff2bed INFO: processing $infile\n";

parse_gff($infile);

feature_summary($fstat,$workdir);

features2bed($feat,$fstat,"xxx",$workdir,$bn,undef);


sub usage {
  print <<EOF;
Convert GFF3 to BED12

usage: $0 [options]
program specific options:                                     default:
 -chr     <string>  specify chromosome name (1st GFF column)  ($chr_name)
 -feature <string>  specify feature to extract                ($feature)
 -i       <string>  input file (GFF3)
 -wd      <string>  working directory (path to GFF file)      ($workdir)
 -help              rint this information
EOF
exit;
}
