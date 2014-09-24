#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-09-24 14:17:55 egg>
#
# Construct UCSC assembly hub
# Display of novel genome sequences on the UCSC Genome Hub.
#
# usage: assembly_hub_constructer.pl -i file.fa -a /output/path -b http://www.test.com/test/
#
# ***********************************************************************
# *  Copyright notice
# *
# *  Copyright 2014 Michael Thomas Wolfinger <michael@wolfinger.eu>
# *                 Florian Michel Eggenhofer <florian@eggenhofer.eu>
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
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use ViennaNGS::UCSC qw( make_assembly_hub );

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
my ($fasta_file_path,$assembly_hub_destination_path,$base_URL,$log_path);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
Getopt::Long::config('no_ignore_case');
&usage() unless GetOptions("i=s"           => \$fasta_file_path,
			   "a=s"           => \$assembly_hub_destination_path,
                           "b=s"           => \$base_URL,
			   "d"             => sub{$debug=1},
                           "-help"         => \&usage,
                           "v");

die "ERROR in [assembly_hub_constructer.pl]: no input fasta file provided" unless(defined $fasta_file_path);
die "ERROR in [assembly_hub_constructer.pl]: input fasta file not found" unless (-e $fasta_file_path);
unless ($infile =~ /^\// || $infile =~ /^\.\//){$infile = "./".$infile;}
unless ($assembly_hub_destination_path =~ /\/$/){$assembly_hub_destination_path.= "/";}
$log_path = $assembly_hub_destination_path . "Log"
#unless (-d $destdir){$cmd = "mkdir -p $destdir"; system($cmd);}
#($basename,$dir,$ext) = fileparse($infile,qr/\.[^.]*/);

make_assembly_hub $fasta_file_path $assembly_hub_destination_path $base_URL $log_path;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#


sub usage {
  print <<EOF;
Calculate normalized expression values (TPM/RPKM) for multicov data

usage: $0 [options]
program specific options:                                default:
 -i      <file>                                          ($infile)
 -a      <path>   output directory                       ($destdir)
 -b      <url>    URL - network location for upload to UCSC                            ($readlength)
 -d               debug output                           ($debug)
 -help            print this information

EOF
exit;
}
