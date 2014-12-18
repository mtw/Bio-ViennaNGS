#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-18 17:14:26 mtw>

use Bio::ViennaNGS::BamStatSummary;
use Data::Dumper;

my @files = ();

#my $bamfile = "/home/mescalin/mtw/Perl/MyModules/Bio-ViennaNGS/share/data1/SRR794848.8k.bam";
my $bamfile = "/scratch/mtw/tmp/J00027/splits.bam/E1_R1.uniq.splits.bam";
push (@files, $bamfile);
my $bamfile = "/scratch/mtw/tmp/encode.bam";
push (@files, $bamfile);
#my $bamfile = "/scratch/mtw/tmp/encode30.bam";
#push (@files, $bamfile);

my $out = "/tmp/";

#$bss1->stat_singleBam();

my $bamsummary = Bio::ViennaNGS::BamStatSummary->new(files => \@files,
						     outpath => $out,
						    );

$bamsummary->populate_data();

$bamsummary->populate_countStat();

$bamsummary->dump_countStat("csv");
#print Dumper(\@files);
#print Dumper($bamsummary);
