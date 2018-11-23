#-*-Perl-*-

use strict;
use warnings;
use Path::Class;
use File::Share ':all';
use Data::Dumper;
use FindBin qw($Bin);
use constant TEST_COUNT => 4;

use lib "$Bin/../lib", "$Bin/../blib/lib", "$Bin/../blib/arch";

BEGIN {
    # include Test.pm from 't' dir in case itis not installed
    eval { require Test::More; };
    if ($@) {
      use lib 't';
    }
    use Test::More tests => TEST_COUNT;
}

use Bio::ViennaNGS::FeatureInterval;

{
  my $chr = "chrX";
  my $start = 99999;
  my $end = 123456;
  my @arg = (chromosome => $chr, start => $start, end => $end);
  my $FI1 = new_ok('Bio::ViennaNGS::FeatureInterval'=> \@arg);
  ok($FI1->_length == 23456, 'length of Featureinterval (base 0)');
  $start = 100;
  $end = 101;
  @arg = (chromosome=>$chr,start=>$start,end=>$end,base=>1);
  my $FI2 = new_ok('Bio::ViennaNGS::FeatureInterval'=> \@arg);
  ok($FI2->_length == 1, 'length of Featureinterval (base 1)');
  
}

