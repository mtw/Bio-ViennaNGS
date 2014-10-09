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
# *                 Florian Eggenhofer <florian.eggenhofer@univie.ac.at>
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
use Bio::ViennaNGS::UCSC qw( make_assembly_hub );

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
my ($fasta_file_path,$assembly_hub_destination_path,$base_URL,$assembly_hub_return_value,$log_path);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
Getopt::Long::config('no_ignore_case');
&usage() unless GetOptions("i=s"           => \$fasta_file_path,
			   "a=s"           => \$assembly_hub_destination_path,
                           "b=s"           => \$base_URL,
                           "-help"         => \&usage,
                           "v");

die "ERROR in [assembly_hub_constructer.pl]: no input fasta file provided" unless(defined $fasta_file_path);
die "ERROR in [assembly_hub_constructer.pl]: input fasta file not found" unless (-e $fasta_file_path);
die "ERROR in [assembly_hub_constructer.pl]: no assembly hub destination path provided" unless(defined $assembly_hub_destination_path);
die "ERROR in [assembly_hub_constructer.pl]: assembly hub destination path not writable" unless (-e $assembly_hub_destination_path);
die "ERROR in [assembly_hub_constructer.pl]: no URL (network location for upload to UCSC) provided" unless(defined $base_URL);
unless ($assembly_hub_destination_path =~ /\/$/){$assembly_hub_destination_path.= "/";}
$log_path = $assembly_hub_destination_path . "Log/";
unless (-d $log_path){my $cmd = "mkdir -p  $log_path"; system($cmd);}

$assembly_hub_return_value = make_assembly_hub($fasta_file_path,$assembly_hub_destination_path,$base_URL,$log_path);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#


sub usage {
  print <<EOF;
Calculate normalized expression values (TPM/RPKM) for multicov data

usage: $0 [options]
program specific options:                                default:
 -i      <file>                                          ($fasta_file_path)
 -a      <path>   output directory                       ($assembly_hub_destination_path)
 -b      <url>    URL - network location for upload to UCSC                            ($base_URL)
 -help            print this information

EOF
exit;
}
