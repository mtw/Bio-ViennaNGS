#-*-Perl-*-

use strict;
use warnings;
use Path::Class;
use File::Share ':all';
use Data::Dumper;
use FindBin qw($Bin);
use constant TEST_COUNT => 1;

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

use Bio::ViennaNGS::MinimalFeature;

{
  my $chr = "chrX";
  my $start = 99999;
  my $end = 123456;
  my $strand = "+";
  my @arg = (chromosome => $chr, start => $start, end => $end, strand => $strand);
  my $MF = new_ok('Bio::ViennaNGS::MinimalFeature'=> \@arg);

}

