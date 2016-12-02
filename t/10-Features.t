#-*-Perl-*-

use strict;
use warnings;
use Path::Class;
use File::Share ':all';
use Data::Dumper;
use FindBin qw($Bin);
use constant TEST_COUNT => 17;

use lib "$Bin/../lib", "$Bin/../blib/lib", "$Bin/../blib/arch";

BEGIN {
    # include Test.pm from 't' dir in case itis not installed
    eval { require Test::More; };
    if ($@) {
      use lib 't';
    }
    use Test::More tests => TEST_COUNT;
    print Dumper($Bin);
}

use Bio::ViennaNGS::Feature;
use Bio::ViennaNGS::FeatureChain;

{
  my $chr = "NC_012345.6";
  my $start1 = 1100; my $start2 = 2345; my $start3 = 2987; my $start4 = 5437;
  my $end1 = 1346; my $end2 = 2544; my $end3 = 3076; my $end4 = 6789;
  my $name1 = "feat1"; my $name2 = "feat2"; my $name3 = "feat3"; my $name4 = "feat4";
  my $strand = "+";
  my @arg1 = (chromosome => $chr, start => $start1, end => $end1, strand => $strand, name => $name1);
  my @arg2 = (chromosome => $chr, start => $start2, end => $end2, strand => $strand, name => $name2);
  my @arg3 = (chromosome => $chr, start => $start3, end => $end3, strand => $strand, name => $name3);
  my @arg4 = (chromosome => $chr, start => $start4, end => $end4, strand => $strand, name => $name4);
  my $Feature1 = new_ok('Bio::ViennaNGS::Feature'=> \@arg1);
  my $Feature2 = new_ok('Bio::ViennaNGS::Feature'=> \@arg2);
  my $Feature3 = new_ok('Bio::ViennaNGS::Feature'=> \@arg3);
  my $Feature4 = new_ok('Bio::ViennaNGS::Feature'=> \@arg4);
  ok($Feature1->_length == 246, 'length of Feature1');
  ok($Feature2->_length == 199, 'length of Feature2');
  ok($Feature3->_length == 89, 'length of Feature3');
  ok($Feature4->_length == 1352, 'length of Feature4');

  my @arg_fc = (type => "test", chain => [$Feature1,$Feature2]);
  my $FC1 = new_ok("Bio::ViennaNGS::FeatureChain" => \@arg_fc);

  ok($FC1->_entries == 2, "elements in FeatureChain->chain ArrayRef");

  # check if the object implements all methods
  can_ok($FC1, qw(count_entries sort_chain_ascending as_bed12_line as_bed6_array));

  # check if $FC->chain is an ArrayRef to 'Bio::ViennaNGS::Bed'
  isa_ok($FC1->chain, 'ARRAY');

  # add a feature to the chain
  $FC1->add($Feature4);
  $FC1->count_entries();
  ok($FC1->_entries == 3, "elements in FeatureChain->chain ArrayRef");

  # and even another one that starts before the last one
  $FC1->add($Feature3);
  $FC1->count_entries();
  ok($FC1->_entries == 4, "elements in FeatureChain->chain ArrayRef");
  $FC1->sort_chain_ascending();

  # check if the chain is sorted by Features' start coordinates
  cmp_ok(${$FC1->chain}[2]->start, '<', ${$FC1->chain}[3]->start, 'chain sorting');

  my $ft = $FC1->pop();
  isa_ok( $ft, 'Bio::ViennaNGS::Feature');
  $FC1->count_entries();
  ok($FC1->_entries == 3, "elements in FeatureChain->chain ArrayRef");
}
