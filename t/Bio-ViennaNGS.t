# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl ViennaNGS.t'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;
use IPC::Cmd qw(can_run);
use File::Share ':all';
use Test::More;
use Data::Dumper;

my $data1_fa  = dist_file('Bio-ViennaNGS','NC_000913.3.30k.fa');
my $data1_gff = dist_file('Bio-ViennaNGS','NC_000913.3.30k.gff');

BEGIN { use_ok('Bio::ViennaNGS') };
BEGIN { use_ok('Bio::ViennaNGS::Fasta') };
BEGIN { use_ok('Bio::ViennaNGS::AnnoC') };
BEGIN { use_ok('Bio::ViennaNGS::Util') };
BEGIN { use_ok('Bio::ViennaNGS::UCSC') };
BEGIN { use_ok('Bio::ViennaNGS::Feature') };
BEGIN { use_ok('Bio::ViennaNGS::MinimalFeature') };
BEGIN { use_ok('Bio::ViennaNGS::FeatureChain') };

ok( defined(can_run('cat')), 'cat not found');
ok( defined(can_run('awk')), 'awk not found');
ok( defined(can_run('bedToBigBed')), 'bedToBigBed not found'); 
ok( defined(can_run('genomeCoverageBed')), 'genomeCoverageBed not found');
ok( defined(can_run('bedGraphToBigWig')), 'bedGraphToBigWig not found');
ok( defined(can_run('faToTwoBit')), 'faToTwoBit not found');
ok( defined(can_run('bedtools')), 'bedtools not found');
ok( defined(can_run('sortBed')), 'sortBed not found');
#########################

my $f = Bio::ViennaNGS::Fasta->new( fa => $data1_fa, );
ok($f);


#my @data1_ids = $f->fastaids;
#ok (scalar @data1_ids, 1);
#my $ps = $f->primaryseq;






# Insert your test code below, the Test::More module is use()ed here so read
# its man page ( perldoc Test::More ) for help writing this test script.

done_testing;
