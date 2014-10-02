#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-10-02 15:01:25 mtw>
#
# Compute normalized expression data in TPM/RPKM from (raw) read
# counts in multicov.
# TPM reference: Wagner et al, Theory Biosci. 131(4), pp 281-85 (2012)
#
# usage: normalize_multicov.pl -i file.multicov -rl 50 -o /output/path
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
use Getopt::Long;
use Data::Dumper;
use File::Basename;
use Bio::ViennaNGS qw( parse_multicov write_multicov featCount_data computeTPM);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
my ($infile,$FC,$FC_sample,$conds,$totreads,$i,$sample,$basename,$dir,$ext,$cmd);
my $readlength    = 100;
my $debug         = 0;
my $destdir       = "./";

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
Getopt::Long::config('no_ignore_case');
&usage() unless GetOptions("i=s"           => \$infile,
			   "rl=s"          => \$readlength,
			   "d"             => sub{$debug=1},
			   "o=s"           => \$destdir,
                           "-help"         => \&usage,
                           "v");

die "ERROR in [normalize_multicov]: no input multicov file provided" unless(defined $infile);
die "ERROR in [normalize_multicov]: input multicov file not found" unless (-e $infile);
unless ($infile =~ /^\// || $infile =~ /^\.\//){$infile = "./".$infile;}
unless ($destdir =~ /\/$/){$destdir .= "/";}
unless (-d $destdir){$cmd = "mkdir -p $destdir"; system($cmd);}
($basename,$dir,$ext) = fileparse($infile,qr/\.[^.]*/);

# parse multicov file; populate @featCount
$conds = parse_multicov($infile);
$FC = featCount_data();
#print Dumper($FC);
#print "parsing $infile: $conds conditions found\n";

# extract hash for a specific sample from $FC
for ($i=0;$i<$conds;$i++){
  my $meanTPM;
  $FC_sample = @$FC[$i];
  $meanTPM = computeTPM($FC_sample, $readlength);
  if ($debug == 1) {
    print "===> sample $i:\n";
    # print Dumper($FC_sample);
  }
  #$totalreads = totalreads(%$FCsample);
  if ($debug == 1){print  "\tmean TPM = $meanTPM\n";}
}

# write multicov file based on TPM data in @featCount
write_multicov("TPM", $destdir, $basename);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#


sub usage {
  print <<EOF;
Calculate normalized expression values (TPM/RPKM) for multicov data

usage: $0 [options]
program specific options:                                default:
 -i      <file>                                          ($infile)
 -o      <path>   output directory                       ($destdir)
 -rl     <int>    read length                            ($readlength)
 -d               debug output                           ($debug)
 -help            print this information

EOF
exit;
}
