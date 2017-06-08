# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl ViennaNGS.t'

#########################

# change 'tests => 3' to 'tests => last_test_to_print';

use strict;
use warnings;
#use IPC::Cmd qw(can_run);
use File::Share ':all';
use Test::More;
use Data::Dumper;

my $data1_fa  = dist_file('Bio-ViennaNGS','NC_000913.3.30k.fa');
my $data1_gff = dist_file('Bio-ViennaNGS','NC_000913.3.30k.gff');

print Dumper($data1_fa);
print Dumper($data1_gff);

BEGIN { use_ok('Bio::ViennaNGS') };
BEGIN { use_ok('Bio::ViennaNGS::AnnoC') };
BEGIN { use_ok('Bio::ViennaNGS::Bam') };
BEGIN { use_ok('Bio::ViennaNGS::Bed') };
BEGIN { use_ok('Bio::ViennaNGS::BedGraphEntry') };
BEGIN { use_ok('Bio::ViennaNGS::Expression') };
BEGIN { use_ok('Bio::ViennaNGS::ExtFeature') };
BEGIN { use_ok('Bio::ViennaNGS::Fasta') };
BEGIN { use_ok('Bio::ViennaNGS::Feature') };
BEGIN { use_ok('Bio::ViennaNGS::FeatureChain') };
BEGIN { use_ok('Bio::ViennaNGS::FeatureLine') };
BEGIN { use_ok('Bio::ViennaNGS::FeatureInterval') };
BEGIN { use_ok('Bio::ViennaNGS::FeatureIO') };
BEGIN { use_ok('Bio::ViennaNGS::MinimalFeature') };
BEGIN { use_ok('Bio::ViennaNGS::Peak') };
BEGIN { use_ok('Bio::ViennaNGS::SpliceJunc') };
BEGIN { use_ok('Bio::ViennaNGS::Tutorial') };
BEGIN { use_ok('Bio::ViennaNGS::UCSC') };
BEGIN { use_ok('Bio::ViennaNGS::Util') };

#ok( defined(can_run('cat')), 'cat not found');
#ok( defined(can_run('awk')), 'awk not found');
#ok( defined(can_run('bedToBigBed')), 'bedToBigBed not found'); 
#ok( defined(can_run('genomeCoverageBed')), 'genomeCoverageBed not found');
#ok( defined(can_run('bedGraphToBigWig')), 'bedGraphToBigWig not found');
#ok( defined(can_run('faToTwoBit')), 'faToTwoBit not found');
#ok( defined(can_run('bedtools')), 'bedtools not found');
#ok( defined(can_run('sortBed')), 'sortBed not found');
#########################

my $chr = "chr1";
my $start = 100;
my $end   = 200;
my $strand = "-";
my $value = 234433.434354;

#my $f = Bio::ViennaNGS::Fasta->new( fa => $data1_fa, );
#ok($f);

my $featureinterval = Bio::ViennaNGS::FeatureInterval->new(
   chromosome =>$chr,
   start => $start,
   end => $end,
); 
ok($featureinterval);

my $mf = Bio::ViennaNGS::MinimalFeature->new(
   chromosome => $chr,
   start => $start,
   end => $end,
   strand => $strand,
); 
ok($mf);

my $bge = Bio::ViennaNGS::BedGraphEntry->new(
   chromosome => $chr,
   start => $start,
   end => $end,
   dataValue => $value,
);
ok($bge);
	


#my @data1_ids = $f->fastaids;
#ok (scalar @data1_ids, 1);
#my $ps = $f->primaryseq;


# Insert your test code below, the Test::More module is use()ed here so read
# its man page ( perldoc Test::More ) for help writing this test script.

done_testing;
