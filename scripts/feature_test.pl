#!/usr/bin/env perl
# Last changed Time-stamp: <2014-11-06 23:11:55 mtw>

use Bio::ViennaNGS::FeatureIO;
use Bio::ViennaNGS::Feature;
use Bio::ViennaNGS::MinimalFeature;
use Bio::ViennaNGS::FeatureLine;
use Bio::ViennaNGS::FeatureChain;
use Bio::ViennaNGS::ContainerFeature;
use Moose;
use MooseX::InstanceTracking;
use Data::Dumper;


my $chr = "chr2";
my $start = 200;
my $end = 400;
my $strand = ".";
my $name1 = "feat1";
my $name2 = "feat2";

my $mf = Bio::ViennaNGS::MinimalFeature->new(chr => $chr,
					     start => $start,
					     end => $end,
					     strand => $strand,
					    );
my $feature1 = Bio::ViennaNGS::Feature->new(
					    chr => $chr,
					    start => $start,
					    end => $end,
					    strand => $strand,
					    name => $name1,
					   );
my $feature2 = Bio::ViennaNGS::Feature->new(
					    chr => $chr,
					    start => "600",
					    end => "800",
					    strand => "+",
					    name => $name2,
					   );

my $fc = Bio::ViennaNGS::FeatureChain->new(
					   type => "test",
					  );
$fc->add($name1 => \$feature1);
$fc->add($name2 => \$feature2);

print Dumper($fc->lookup($name2));


#print Dumper($feature1);
#print Dumper($feature2);
print Dumper($fc);
print "===================================\n";
#print "instances of Bio::ViennaNGS::FeatureChain:\n";
#my $xxx = Bio::ViennaNGS::FeatureChain->meta->instances;
#print Dumper($xxx);

my @feature_instances = Bio::ViennaNGS::Feature->meta->instances;
foreach my $f (@feature_instances){
  print "++++++++++++++++++\n";
  print Dumper ($f);
 print "\n\n";
}

my $cf1 = Bio::ViennaNGS::FeatureIO->new();

my $cf2 = Bio::ViennaNGS::ContainerFeature->new(
					       chr => $feature1->chr,
					       start => $feature1->start,
					       end => $feature1->end,
					       strand => $feature1->strand,
					       parent => $cf1, 
					   
					      );
print Dumper($cf2);

