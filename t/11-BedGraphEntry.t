#-*-Perl-*-

use strict;
use warnings;
use Path::Class;
use File::Share ':all';
use Data::Dumper;
use FindBin qw($Bin);
use constant TEST_COUNT => 2;

use lib "$Bin/../lib", "$Bin/../blib/lib", "$Bin/../blib/arch";

BEGIN {
    # include Test.pm from 't' dir in case itis not installed
    eval { require Test::More; };
    if ($@) {
      use lib 't';
    }
    use Test::More tests => TEST_COUNT;
}

use Bio::ViennaNGS::BedGraphEntry;

{
  my $chr = "chrX";
  my $start = 99999;
  my $end = 123456;
  my $value = 234433.434354;
  my @arg = (chromosome => $chr, start => $start, end => $end, dataValue => $value);
  my $BGE = new_ok('Bio::ViennaNGS::BedGraphEntry'=> \@arg);
  ok($BGE->_length == 23456, 'length of underlying Featureinterval (base 0)');
}

