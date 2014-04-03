#!/usr/bin/env perl -w
# -*-CPerl-*-
# Last changed Time-stamp: <2013-09-04 17:02:56 mtw>
#
# Convert GFF3 to BED12; produce separate BED files for each gbkey (CDS/tRNA/etc)
#
# usage: gff2bed.pl -i input.gff
#
# ***********************************************************************
# *  Copyright notice
# *
# *  Copyright 2013 Michael Thomas Wolfinger <michael@wolfinger.eu>
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
use local::lib;
use ViennaNGS::AnnoC;

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
make_beds_from_features($feat,$fstat);


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
sub make_beds_from_features {
  my ($f,$stat) = @_;
  my ($chrom, $chrom_start, $chrom_end, $name, $score, $strand, $thick_start, $thick_end);
  my ($reserved, $block_count, $block_sizes, $block_starts);
  
 # print "in make_beds:...\n";
 # print Dumper (\$f);
  # 1. loop over all feature types found in GFF3
  foreach my $gbkey (keys %$stat) {
    next if ($gbkey eq 'total' || $gbkey eq 'Src' || $gbkey eq 'accession' || $gbkey eq 'origin');
    my $bednameuns = $workdir.$bn.".".$gbkey.".uns.bed";
    my $bedname    = $workdir.$bn.".".$gbkey.".bed";
    print "Processing gbkey $gbkey...into $bedname\n";
    open (BEDFILE, "> $bednameuns") or die "ERROR: cannot open $bednameuns for writing\n";
    # 2. print unsorted info from DS extracted from GFF3 as BED
    foreach my $uid (keys %$f){
      next unless ($$f{$uid}->{gbkey} eq $gbkey);
      my @bedline = ();
      $chrom        = $$f{$uid}->{seqid};
      $chrom_start  = $$f{$uid}->{start};
      $chrom_start--; # BED is 0-based
      $chrom_end    = $$f{$uid}->{end};
      $name         = $$f{$uid}->{name};
      $score        = $$f{$uid}->{score};
      $strand       = $$f{$uid}->{strand} == -1 ? '-' : '+'; #default to +
      $thick_start  = $chrom_start;
      $thick_end    = $chrom_end;
      $reserved     = 0; 
      $block_count  = 1;
      $block_sizes  = $$f{$uid}->{length}.",";
      $block_starts = "0,";
      @bedline = join ("\t", ($chrom,$chrom_start,$chrom_end,$name,$score,$strand,$thick_start,$thick_end,$reserved,$block_count,$block_sizes, $block_starts));
      print BEDFILE "@bedline\n";
    }
    close BEDFILE;

    my $cmdl = $sortBed." -i ".$bednameuns." > ".$bedname;
    system($cmdl);
    unlink($bednameuns);
  }
}

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
