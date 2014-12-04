#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-04 22:20:43 mtw>
#
# Create a new genome database for the UCSC genome browser. Based on
# http://genomewiki.ucsc.edu/index.php/Building_a_new_genome_database
#
# usage: newUCSCdb.pl -db <DBNAME> -fa <FASTA> -genome <GENOME> -tax <TAXID>
#                     -sname <SCIENTIFIC_NAME> -assembly <ASSEMBLY> -clade <CLADE>
#                     -pos <DEFAULT_POSITION> -source <SOURCE> -key <INT>
#
#
# EXAMPLE 1:
# newUCSCdb.pl -db eColiMC4100 -fa NC_012759.fna -tax 595496 -genome
# "Escherichia coli MC4100" -assembly "2013-09-09" -pos
# "NC_012759.1:230000-867000" -key 204 -sname "Escherichia coli BW2952
# uid 59391" -source "based on RefSeq NC_012759.1" -clade "bacteria"
# -h
#
# EXAMPLE 2:
# newUCSCdb.pl -db stenMaltK279a -tax 522373 -genome "Stenotrophomonas
# maltophilia K279a" -assembly "2008-06-10" -pos "NC_010943.1:
# 100000-200000" -fa NC_010943.fna -sname "Stenotrophomonas
# maltophilia K279a uid 61647" -key 223 -source "based on RefSeq
# NC_010943.1" -clade "bacteria" -h

# ***********************************************************************
# *  Copyright notice
# *
# *  Copyright 2013 Michael T. Wolfinger <michael@wolfinger.eu>
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
use POSIX qw(strftime);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my $timestamp    = strftime("%Y%m%d-%H%M", localtime(time));       # timestamp
my $awk          = "awk";                        # path to awk
my $hgsql        = "hgsql";                      # path to hgsql
my $hgFakeAgp    = "hgFakeAgp ";                 # path to hgFakeAgp
my $faToTwoBit   = "faToTwoBit ";                # path to faToTwoBit
my $twoBitInfo   = "twoBitInfo ";                # path to twoBitInfo
my $hgLoadSqlTab = "hgLoadSqlTab ";              # path to hgLoadSqlTab
my $hgGcPercent  = "hgGcPercent ";
my $hgLoadWiggle = "hgLoadWiggle ";
my $hgGoldGapGl  = "hgGoldGapGl ";
my $log          = "newUCSCdb.log";
my $db = my $fa = my $taxid = my $genome = my $sname = my $clade = undef;
my $assembly = my $position = my $orderKey = my $source = undef;


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

Getopt::Long::config('no_ignore_case');
&usage() unless GetOptions("db=s"       => \$db,
                           "fa=s"       => \$fa,
                           "tax=s"      => \$taxid,
			   "genome=s"   => \$genome,
			   "sname=s"    => \$sname,
			   "assembly=s" => \$assembly,
			   "pos=s"      => \$position,
			   "key=i"      => \$orderKey,
			   "source=s"   => \$source,
			   "clade=s"    => \$clade,
                           "-help"      => \&usage,
                           "v");

unless ((defined $db) && (defined $fa) && (defined $taxid) && (defined $genome)) {
  die "Please provide command-line arguments for -db / -fa / -tax / -genome \n";
}
unless (defined $sname) {
  $sname = $genome;
}
unless (defined $assembly) {
  $assembly = "Jan 2041";
}
unless (defined $position) {
  $position = "chr1:500000-1000000";
}
unless (defined $orderKey) {
  $orderKey = 9999;
}
unless (defined $source) {
  $source = "new genome assembly";
}
unless (defined $clade) {
  $clade = "other"; # mammal / vertebrate / plants / bacteria
}

make_ucsc_db();


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

sub make_ucsc_db {
  my $cmd;
  my @log = ();
  open my $logfile, ">", $log;

  # log header
  push @log, "***   START $timestamp\n";

  # create agp file
  $cmd = $hgFakeAgp." -minContigGap=1 ".$fa." $db.agp\n";
  system($cmd);
  push @log, $cmd; push @log, "\n";

  # convert fasta to 2bit
  $cmd = "$faToTwoBit $fa $db.2bit\n";
  $cmd .= "mkdir -p /gbdb/$db/html\n";
  $cmd .= "cp $db.2bit /gbdb/$db/\n";
  system($cmd);
  push @log, $cmd; push @log, "\n";

  # create chromInfo table
  $cmd = "$twoBitInfo $db.2bit stdout | sort -k2nr > chrom.sizes\n";
  $cmd .= "mkdir -p bed/chromInfo\n";
  $cmd .= "awk \' {printf \"%s\\t%d\\t/gbdb/$db/$db.2bit\\n\", \$1, \$2}\' chrom.sizes > bed/chromInfo/chromInfo.tab\n";
  system($cmd);
  push @log, $cmd; push @log, "\n";

  # start a new database
  $cmd = "hgsql -e \"CREATE DATABASE IF NOT EXISTS $db;\" mysql\n";
  system($cmd);
  push @log, $cmd; push @log, "\n";

  # load grp table;
  $cmd = "hgsql $db < \$HOME/kent/src/hg/lib/grp.sql\n";
  system($cmd);
  push @log, $cmd; push @log, "\n";

  # load chromInfo table;
  $cmd = "$hgLoadSqlTab $db chromInfo \$HOME/kent/src/hg/lib/chromInfo.sql bed/chromInfo/chromInfo.tab\n";
  system($cmd);
  push @log, $cmd; push @log, "\n";

  # load gold and gap tables
  $cmd = "$hgGoldGapGl $db $db.agp\n";
  system($cmd);
  push @log, $cmd; push @log, "\n";

  # generate gc5Base data and load the table
  $cmd = "mkdir -p bed/gc5Base\n";
  $cmd .= "$hgGcPercent -wigOut -doGaps -file=gc5.wig -win=5 -verbose=0 $db $db.2bit \n";
  $cmd .= "wigEncode  gc5.wig bed/gc5Base/gc5Base.wig bed/gc5Base/gc5Base.wib \n";
  $cmd .= "$hgLoadWiggle $db gc5Base bed/gc5Base/gc5Base.wig\n";
  $cmd .= "mkdir -p /gbdb/$db/wib\n";
  $cmd .= "cp bed/gc5Base/gc5Base.wib /gbdb/$db/wib\n";
  system($cmd);
  push @log, $cmd; push @log, "\n";

  # create dbDbInsert.sql
  open my $dbdb, ">", "dbDbInsert.sql";
  print $dbdb "INSERT INTO dbDb\n";
  print $dbdb "  (name, description, nibPath, organism,\n";
  print $dbdb "  defaultPos, active, orderKey, genome, scientificName,\n";
  print $dbdb "  htmlPath, hgNearOk, hgPbOk, sourceName, taxId)\n";
  print $dbdb "VALUES\n";
  print $dbdb "  (\"$db\", \"$assembly\", \"/gbdb/$db\", \"$genome\",\n";
  print $dbdb "  \"$position\", 1, $orderKey, \"$genome\", \"$sname\",\n";
  print $dbdb "  \"/gbdb/$db/html/description.html\", 0, 0, \"$source\", $taxid);\n";
  close $dbdb;

  # load dbDbInsert.sql into database
  $cmd = "hgsql hgcentral < dbDbInsert.sql\n";
  system($cmd);
  push @log, $cmd; push @log, "\n";

  # create defaultDb entry
  $cmd = "hgsql hgcentral -e \'INSERT INTO defaultDb (genome, name) VALUES (\"$genome\",\"$db\");\'\n";
  system($cmd);
  push @log, $cmd; push @log, "\n";

  # create genomeClade entry
  $cmd = "hgsql hgcentral -e 'INSERT INTO genomeClade (genome, clade, priority) VALUES (\"$genome\", \"$clade\", $orderKey);\'\n";
  system($cmd);
  push @log, $cmd; push @log, "\n";

  # log footer
  push @log, "*** END $timestamp\n";

  print $logfile @log;
  close $log;

}


sub usage {
  print <<EOF;
Depoly a new UCSC genome browser database

usage: $0 [options]
program specific options:                                    default:
 -db       <string>  new database name                        ($db)
 -fa       <string>  (multi) fasta file holding sequence      ($fa)
 -tax      <string>  taxonomy ID                              ($taxid)
 -genome   <string>  genome name                              ($genome)
 -sname    <string>  scientific name                          ($sname)
 -assembly <string>  assembly                                 ($assembly)
 -pos      <string>  default position                         ($position)
 -key      <int>     order key                                ($orderKey)
 -source   <string>  genome/assembly/annotation source        ($source)
 -clade    <string>  genome clade (mammal,plants,bateria,...) ($clade)
 -help               print this information
EOF
exit;
}
