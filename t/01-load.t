#-*-Perl-*-

use strict;
use warnings;
use Test::More tests => 22;

BEGIN {
      use_ok('Bio::ViennaNGS') ||  print "Bail out! Cannot load Bio::ViennaNGS\n";
      use_ok('Bio::ViennaNGS::AnnoC')  ||  print "Bail out! Cannot load Bio::ViennaNGS::AnnoC\n";	
      use_ok('Bio::ViennaNGS::Bam')  ||  print "Bail out! Cannot load Bio::ViennaNGS::Bam\n";
      use_ok('Bio::ViennaNGS::Bed')  ||  print "Bail out! Cannot load Bio::ViennaNGS::Bed\n";
      use_ok('Bio::ViennaNGS::BedGraphEntry')  ||  print "Bail out! Cannot load Bio::ViennaNGS::BedGraphEntry\n";
      use_ok('Bio::ViennaNGS::Expression')  ||  print "Bail out! Cannot load Bio::ViennaNGS::Expression\n";
      use_ok('Bio::ViennaNGS::ExtFeature')  ||  print "Bail out! Cannot load Bio::ViennaNGS::ExtFeature\n";
      use_ok('Bio::ViennaNGS::Fasta')  ||  print "Bail out! Cannot load Bio::ViennaNGS::Fasta\n";
      use_ok('Bio::ViennaNGS::Feature')  ||  print "Bail out! Cannot load Bio::ViennaNGS::Feature\n";
      use_ok('Bio::ViennaNGS::FeatureBase')  ||  print "Bail out! Cannot load Bio::ViennaNGS::FeatureBase\n";
      use_ok('Bio::ViennaNGS::FeatureChain')  ||  print "Bail out! Cannot load Bio::ViennaNGS::FeatureChain\n";
      use_ok('Bio::ViennaNGS::FeatureLine')  ||  print "Bail out! Cannot load Bio::ViennaNGS::FeatureLine\n";
      use_ok('Bio::ViennaNGS::FeatureInterval')  ||  print "Bail out! Cannot load Bio::ViennaNGS::Featureinterval\n";
      use_ok('Bio::ViennaNGS::FeatureIntervalN')  ||  print "Bail out! Cannot load Bio::ViennaNGS::FeatureintervalN\n";
      use_ok('Bio::ViennaNGS::FeatureIO')  ||  print "Bail out! Cannot load Bio::ViennaNGS::FeatureIO\n";
      use_ok('Bio::ViennaNGS::MinimalFeature')  ||  print "Bail out! Cannot load Bio::ViennaNGS::MinimalFeature\n";
      use_ok('Bio::ViennaNGS::Peak')  ||  print "Bail out! Cannot load Bio::ViennaNGS::Peak\n";
      use_ok('Bio::ViennaNGS::SpliceJunc')  ||  print "Bail out! Cannot load Bio::ViennaNGS::SpliceJunc\n";
      use_ok('Bio::ViennaNGS::Subtypes')  ||  print "Bail out! Cannot load Bio::ViennaNGS::Subtypes\n";
      use_ok('Bio::ViennaNGS::Tutorial')  ||  print "Bail out! Cannot load Bio::ViennaNGS::tutorial\n";
      use_ok('Bio::ViennaNGS::UCSC')  ||  print "Bail out! Cannot load Bio::ViennaNGS::UCSC\n";
      use_ok('Bio::ViennaNGS::Util')  ||  print "Bail out! Cannot load Bio::ViennaNGS::Util\n";
}


diag( "Test Bio::ViennaNGS $Bio::ViennaNGS::VERSION, Perl $], $^X" );
diag( "Test Bio::ViennaNGS::AnnoC $Bio::ViennaNGS::AnnoC::VERSION, Perl $], $^X" );
diag( "Test Bio::ViennaNGS::Bam $Bio::ViennaNGS::Bam::VERSION, Perl $], $^X" );
diag( "Test Bio::ViennaNGS::Bed $Bio::ViennaNGS::Bed::VERSION, Perl $], $^X" );
diag( "Test Bio::ViennaNGS::BedGraphEntry $Bio::ViennaNGS::BedGraphEntry::VERSION, Perl $], $^X" );
diag( "Test Bio::ViennaNGS::Expression $Bio::ViennaNGS::Expression::VERSION, Perl $], $^X" );
diag( "Test Bio::ViennaNGS::ExtFeature $Bio::ViennaNGS::ExtFeature::VERSION, Perl $], $^X" );
diag( "Test Bio::ViennaNGS::Fasta $Bio::ViennaNGS::Fasta::VERSION, Perl $], $^X" );
diag( "Test Bio::ViennaNGS::Feature $Bio::ViennaNGS::Feature::VERSION, Perl $], $^X" );
diag( "Test Bio::ViennaNGS::FeatureBase $Bio::ViennaNGS::FeatureBase::VERSION, Perl $], $^X" );
diag( "Test Bio::ViennaNGS::FeatureChain $Bio::ViennaNGS::FeatureChain::VERSION, Perl $], $^X" );
diag( "Test Bio::ViennaNGS::FeatureLine $Bio::ViennaNGS::FeatureLine::VERSION, Perl $], $^X" );
diag( "Test Bio::ViennaNGS::FeatureInterval $Bio::ViennaNGS::FeatureInterval::VERSION, Perl $], $^X" );
diag( "Test Bio::ViennaNGS::FeatureIntervalN $Bio::ViennaNGS::FeatureIntervalN::VERSION, Perl $], $^X" );
diag( "Test Bio::ViennaNGS::FeatureIO $Bio::ViennaNGS::FeatureIO::VERSION, Perl $], $^X" );
diag( "Test Bio::ViennaNGS::MinimalFeature $Bio::ViennaNGS::MinimalFeature::VERSION, Perl $], $^X" );
diag( "Test Bio::ViennaNGS::Peak $Bio::ViennaNGS::Peak::VERSION, Perl $], $^X" );
diag( "Test Bio::ViennaNGS::SpliceJunc $Bio::ViennaNGS::SpliceJunc::VERSION, Perl $], $^X" );
diag( "Test Bio::ViennaNGS::Subtypes $Bio::ViennaNGS::Subtypes::VERSION, Perl $], $^X" );
diag( "Test Bio::ViennaNGS::Tutorial $Bio::ViennaNGS::Tutorial::VERSION, Perl $], $^X" );
diag( "Test Bio::ViennaNGS::UCSC $Bio::ViennaNGS::UCSC::VERSION, Perl $], $^X" );
diag( "Test Bio::ViennaNGS::Util $Bio::ViennaNGS::Util::VERSION, Perl $], $^X" );
