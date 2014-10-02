# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl ViennaNGS.t'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;
use IPC::Cmd qw(can_run);

use Test::More tests => 8;
BEGIN { use_ok('Bio::ViennaNGS') };

ok( defined(can_run('cat')), 'cat not found');
ok( defined(can_run('awk')), 'awk not found');
ok( defined(can_run('bedToBigBed')), 'bedToBigBed not found'); 
ok( defined(can_run('genomeCoverageBed')), 'genomeCoverageBed not found');
ok( defined(can_run('bedGraphToBigWig')), 'bedGraphToBigWig not found');
ok( defined(can_run('bedtools')), 'bedtools not found');
ok( defined(can_run('sortBed')), 'sortBed not found')
#########################

# Insert your test code below, the Test::More module is use()ed here so read
# its man page ( perldoc Test::More ) for help writing this test script.

