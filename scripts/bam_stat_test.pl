#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2014-12-17 00:29:23 mtw>

use Bio::ViennaNGS::BamStatSingle;
use Data::Dumper;

my $bamfile = "/path/to/file.bam";
my %hash = ( 'foo' => { 'aaa' => 'franz',
			'bbb' => 'fritz'},
	     'bar' => { 'aaa' => 'mitzi',
			'bbb' => 'zenzi'}
	   );

my $bss1 = Bio::ViennaNGS::BamStatSingle->new(
					      bam => $bamfile,
					      data => \%hash,
					     );



print Dumper($bss1);
