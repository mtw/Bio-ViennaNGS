#-*-Perl-*-

use strict;
use warnings;
use Path::Class;
use File::Share ':all';
use Data::Dumper;
use FindBin qw($Bin);
use constant TEST_COUNT => 104;

use lib "$Bin/../lib", "$Bin/../blib/lib", "$Bin/../blib/arch";

BEGIN {
    # include Test.pm from 't' dir in case itis not installed
    eval { require Test::More; };
    if ($@) {
      use lib 't';
    }
    use Test::More tests => TEST_COUNT;
}

use Bio::ViennaNGS::FeatureIO;

{
  my $infile_path = dist_dir('Bio-ViennaNGS');
  ok($infile_path);
  my $infile_bed12  = dist_file('Bio-ViennaNGS','data1/NC_000913.CDS.bed');
  ok($infile_bed12);
  my ($infile_bed6) = dist_file('Bio-ViennaNGS','data1/NC_000913.CDS.bed6');
  ok($infile_bed6);

  my @arg1 = (file => $infile_bed12,filetype => 'Bed12',instanceOf => 'Bed');
  my $FIO1 = new_ok("Bio::ViennaNGS::FeatureIO" => \@arg1);

  ok($FIO1->_entries == 30, "elements in ArrayRef");

  # check if the object implements all methods
  can_ok($FIO1, qw(count_entries parse_bedgraph_file parse_bed6_file parse_bed12_file));

  # check if $FIO->data is an ArrayRef to 'Bio::ViennaNGS::Bed'
  isa_ok($FIO1->data, 'ARRAY');

  # check if @{$FIO->data} contains correct objects
  foreach my $bo (@{$FIO1->data}){
    isa_ok( $bo, 'Bio::ViennaNGS::Bed');
  }

  #------
  # test parsing of Bed6 into Bio::ViennaNGS::Feature objects
  my @arg2 = (file => $infile_bed6,filetype => 'Bed6',instanceOf => 'Feature');
  my $FIO2 = new_ok("Bio::ViennaNGS::FeatureIO" => \@arg2);

  ok($FIO2->_entries == 30, "elements in ArrayRef");
  foreach my $bo (@{$FIO2->data}){
    isa_ok( $bo, 'Bio::ViennaNGS::Feature');
  }

  #------
  # test parsing of Bed6 into Bio::ViennaNGS::FeatureChain objects with 1 entry
  my @arg3 =  (file => $infile_bed6,filetype => 'Bed6',instanceOf => 'FeatureChain');
  my $FIO3 = new_ok("Bio::ViennaNGS::FeatureIO" => \@arg3);

  ok($FIO3->_entries == 30, "elements in ArrayRef");
  foreach my $bo (@{$FIO3->data}){
    isa_ok( $bo, 'Bio::ViennaNGS::FeatureChain');
  }

  #------
  # test parsing of Bed6 into Bio::ViennaNGS::FeatureChain object with many entries
  my @arg4 =  (file => $infile_bed6,filetype => 'Bed6',instanceOf => 'FeatureChainBlock');
  my $FIO4 = new_ok("Bio::ViennaNGS::FeatureIO" => \@arg4);

  ok($FIO4->_entries == 1, "elements in ArrayRef");
  foreach my $bo (@{$FIO4->data}){
    isa_ok( $bo, 'Bio::ViennaNGS::FeatureChain');
  }

}

