#-*-Perl-*-

use strict;
use warnings;
use Path::Class;
use File::Share ':all';
use Data::Dumper;
use FindBin qw($Bin);
use constant TEST_COUNT => 204;

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
  my $infile_bg = dist_file('Bio-ViennaNGS','data1/example.bg');
  ok($infile_bg);


  #------
  # test parsing of Bed12 into an array of Bio::ViennaNGS::Bed objects
  my @arg1 = (file => $infile_bed12,filetype => 'Bed12',instanceOf => 'Bed');
  my $FIO1 = new_ok("Bio::ViennaNGS::FeatureIO" => \@arg1);

  # check if the object implements all methods
  can_ok($FIO1, qw(count_entries parse_bedgraph_file parse_bed6_file parse_bed12_file));

   ok($FIO1->_entries == 30, "elements in ArrayRef");
  isa_ok($FIO1->data, 'ARRAY');
  foreach my $bo (@{$FIO1->data}){
    isa_ok( $bo, 'Bio::ViennaNGS::Bed');
  }

  #------
  # test parsing of Bed6 into an array of Bio::ViennaNGS::Feature objects
  my @arg2 = (file => $infile_bed6,filetype => 'Bed6',instanceOf => 'Feature', base => 0);
  my $FIO2 = new_ok("Bio::ViennaNGS::FeatureIO" => \@arg2);

  ok($FIO2->_entries == 30, "elements in ArrayRef");
  foreach my $bo (@{$FIO2->data}){
    isa_ok( $bo, 'Bio::ViennaNGS::Feature');
    ok($bo->base=="0", "Feature element is 0-based");
  }

  #------
  # test parsing of Bed6 into an array of Bio::ViennaNGS::FeatureChain
  # objects, with 1 element each
  my @arg3 =  (file => $infile_bed6,filetype => 'Bed6',instanceOf => 'FeatureChain', base => 0);
  my $FIO3 = new_ok("Bio::ViennaNGS::FeatureIO" => \@arg3);

  ok($FIO3->_entries == 30, "elements in ArrayRef");
  foreach my $bo (@{$FIO3->data}){
    isa_ok( $bo, 'Bio::ViennaNGS::FeatureChain');
    ok($bo->base=="0", "FeatureChain is 0-based");
  }

  #------
  # test parsing of Bed6 into Bio::ViennaNGS::FeatureChain object with many entries
  my @arg4 =  (file => $infile_bed6,filetype => 'Bed6',instanceOf => 'FeatureChainBlock', base => 0);
  my $FIO4 = new_ok("Bio::ViennaNGS::FeatureIO" => \@arg4);

  ok($FIO4->_entries == 1, "elements in ArrayRef");
  foreach my $bo (@{$FIO4->data}){
    isa_ok( $bo, 'Bio::ViennaNGS::FeatureChain');
    while (my $fo = $bo->pop()){
     ok($fo->base=="0", "FeatureChain entry is 0-based");
    }
  }

  #------
  # test parsing of BedGraph into an array of Bio::ViennaNGS::BedGraphEntry objects
  my @arg_bg = (file => $infile_bg, filetype => 'BedGraph', instanceOf => 'BedGraph');
  my $FIO_BG = new_ok("Bio::ViennaNGS::FeatureIO" => \@arg_bg);

  ok(scalar(@{$FIO_BG->data}) == 3, "elements in ArrayRef");
  foreach my $bge (@{$FIO_BG->data}){
    isa_ok( $bge, 'Bio::ViennaNGS::BedGraphEntry');
  }
  my $first = ${$FIO_BG->data}[0];
  ok($first->chromosome eq 'chr1', 'chromosome parsed');
  ok($first->start == 0, 'start parsed');
  ok($first->end == 10, 'end parsed');
  ok($first->dataValue == 1.5, 'value parsed');
}

