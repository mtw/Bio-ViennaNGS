#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-17 17:52:17 mtw>

use Bio::ViennaNGS::BamStatSingle;
use Data::Dumper;
use File::Share ':all';


#my $bamfile = "/home/mescalin/mtw/Perl/MyModules/Bio-ViennaNGS/share/data1/SRR794848.8k.bam";
my $bamfile = "/scratch/mtw/tmp/J00027/splits.bam/E1_R1.uniq.splits.bam";
#my $bamfile = "/scratch/mtw/tmp/encode.bam";
my $bamfile = "/scratch/mtw/tmp/encode30.bam";

#my $bamfile = dist_file('Bio-ViennaNGS','SRR794848.8k.bam');

my %hash = ( 'foo' => { 'aaa' => 'franz',
			'bbb' => 'fritz'},
	     'bar' => { 'aaa' => 'mitzi',
			'bbb' => 'zenzi'}
	   );

my $bss1 = Bio::ViennaNGS::BamStatSingle->new(
					      bam => $bamfile,
					      control_clip => 1,
					     );

$bss1->stat_singleBam();



print Dumper($bss1);
