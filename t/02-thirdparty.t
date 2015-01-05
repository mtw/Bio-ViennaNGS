#!/usr/bin/env perl
use strict;
use warnings;
use Test::More tests => 9;
use IPC::Cmd qw(can_run);

SKIP: {
     skip 'cat seemingly not available',1 
       unless defined(can_run('cat')) ;
     ok( defined(can_run('cat')), 'cat available and executable');
}

SKIP: {
     skip 'awk seemingly not available',1 
       unless defined(can_run('awk')) ;
     ok( defined(can_run('awk')), 'awk available and executable');
}

SKIP: {
     skip 'bedToBigBed seemingly not available',1 
       unless defined(can_run('bedToBigBed')) ;
     ok( defined(can_run('bedToBigBed')), 'bedToBigBed available and executable');
}

SKIP: {
     skip 'genomeCoverageBed seemingly not available',1 
       unless defined(can_run('genomeCoverageBed')) ;
     ok( defined(can_run('genomeCoverageBed')), 'genomeCoverageBed available and executable');
}

SKIP: {
     skip 'bedGraphToBigWig seemingly not available',1 
       unless defined(can_run('bedGraphToBigWig')) ;
     ok( defined(can_run('bedGraphToBigWig')), 'bedGraphToBigWig available and executable');
}

SKIP: {
     skip ('faToTwoBit seemingly is not available',1)
       unless defined(can_run('faToTwoBit'));

     ok( defined(can_run('faToTwoBit')), 'faToTwoBit available and executable') or
       diag ("please install the 'faToTwoBit' utility"); 
}

SKIP: {
     skip ('bedtools seemingly is not available',1)
       unless defined(can_run('bedtools'));
     ok( defined(can_run('bedtools')), 'bedtools available and executable');
}

SKIP: {
     skip ('sortBed seemingly is not available',1)
       unless defined(can_run('sortBed'));
     ok( defined(can_run('sortBed')), 'sortBed available and executable');
}

SKIP: {
     skip ('R seemingly is not available',1)
       unless defined(can_run('R'));
     ok( defined(can_run('R')), 'R available and executable');
}



