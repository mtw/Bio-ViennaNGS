use strict;
use warnings;
use Test::More tests => 7;
use IPC::Cmd qw(can_run);


ok( defined(can_run('cat')), 'cat not found');
ok( defined(can_run('awk')), 'awk not found');
ok( defined(can_run('bedToBigBed')), 'bedToBigBed not found'); 
ok( defined(can_run('genomeCoverageBed')), 'genomeCoverageBed not found');
ok( defined(can_run('bedGraphToBigWig')), 'bedGraphToBigWig not found');

SKIP {
     skip 'faToTwoBit is not available',1 unless {
     defined(can_run('faToTwoBit')) };
    #  ok( defined(can_run('faToTwoBit')), 'faToTwoBit not found');
}
ok( defined(can_run('bedtools')), 'bedtools not found');
ok( defined(can_run('sortBed')), 'sortBed not found');
