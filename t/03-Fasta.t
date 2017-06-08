#-*-Perl-*-

use strict;
use warnings;
use File::Share ':all';
use Data::Dumper;
use FindBin qw($Bin);
use constant TEST_COUNT => 8;

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

use Bio::ViennaNGS::Fasta;

{
  my $data1_fa  = dist_file('Bio-ViennaNGS','data3/all_DENVG_5UTR.fa');
  my @arg1 = (fasta => $data1_fa);
  my $fo1 = new_ok('Bio::ViennaNGS::Fasta'=> \@arg1);
  ok(scalar(@{$fo1->fastaids}) == 5, '# of Fasta entries');

  # check if the object implements all methods
  can_ok($fo1, qw(stranded_subsequence has_sequid));

  # check if $fo->fastaids is an ArrayRef
  isa_ok($fo1->fastaids, 'ARRAY');

  # check if $fo->primaryseqH is a HashRef
  isa_ok($fo1->primaryseqH, 'HASH');

  # check if value of $fo->primaryseqH is a Bio::PrimarySeq::Fasta object
  isa_ok( ${$fo1->primaryseqH}{'NC_001477.1'}, 'Bio::PrimarySeq::Fasta');

  # check if extracted subsequence is correct [+] strand
  ok($fo1->stranded_subsequence("NC_001477.1",17,45,1) eq "GGACCGACAAGAACAGTTTCGAATCGGAA", "check subsequence [+] strand");

   # check if extracted subsequence is correct [-] strand
  ok($fo1->stranded_subsequence("NC_001477.1",71,85,'-') eq "GCTCTCTAATAAAAA", "check subsequence [-] strand");

}
