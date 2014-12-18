#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-18 14:18:02 mtw>

use Bio::ViennaNGS::BamStat;
use Data::Dumper;
use File::Share ':all';

my @files = ();

#my $bamfile = "/home/mescalin/mtw/Perl/MyModules/Bio-ViennaNGS/share/data1/SRR794848.8k.bam";
my $bamfile = "/scratch/mtw/tmp/J00027/splits.bam/E1_R1.uniq.splits.bam";
my $bamfile = "/scratch/mtw/tmp/encode.bam";
#my $bamfile = "/scratch/mtw/tmp/encode30.bam";

#my $bamfile = dist_file('Bio-ViennaNGS','SRR794848.8k.bam');

push (@files, $bamfile);

my $bss1 = Bio::ViennaNGS::BamStatSingle->new(
					      bam => $bamfile,
					      control_clip => 1,
					     );

$bss1->stat_singleBam();



print Dumper($bss1);
print Dumper(\@files);
