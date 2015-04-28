use strict;
use Bio::ViennaNGS::UCSC qw ( make_assembly_hub
                              parse_fasta_header
                              valid_ncbi_accession
                              modify_fasta_header
                              make_group
                              retrieve_chromosome_size
                              write_chromosome_size_file
                              convert_tracks );

use File::Share ':all';
use Test::Files;
use Test::More tests => 8;
use Test::File::Contents;
use Test::Deep;
use Data::Dumper;

my $fasta  = dist_file('Bio-ViennaNGS','data1/NC_000913.3.30k.fa');
print Dumper($fasta);
# Make assembly hub
make_assembly_hub ( $fasta,
                    ".",
                    ".",
                    "http://tbi.univie.ac.at/~egg/assemblyHub_test",
                    "-",
                    "log_asemblyHub.txt" );


# Check if functions work correctly
subtest 'parse_fasta_header' => sub {
  plan tests => 1;
  my @arguments = ("NC_000913.3.30k.fa", "Bio::ViennaNGS::UCSC::make_assembly_hub");
  my $expected = bag ('NC_000913.3', 'scientific name not set');

  cmp_deeply ( parse_fasta_header (@arguments), $expected, 'parse fasta header');
};

subtest 'valid_ncbi_accession' => sub {
  plan tests => 1;
  my $argument = "NC_000913.3";
  my $expected = "NC_000913.3";
is ( valid_ncbi_accession ($argument), $expected, 'valid ncbi accession');
};

subtest 'modify_fasta_header' => sub {
  plan tests => 1;
  my @arguments = ( "NC_000913.3.30k.fa", "assemblyHub/NC_000913.3.fa", "NC_000913.3" );
  my $expected = "1";
is ( modify_fasta_header (@arguments), $expected, 'modify fasta header');
};

subtest 'make_group' => sub {
  plan tests => 1;
  my @arguments = ( "annotation", "Annotation", "1", "0" );
  my $expected =
"name annotation
label Annotation
priority 1
defaultIsClosed 0
";
is ( make_group (@arguments), $expected, 'make group');
};

subtest 'retrieve_chromosome_size' => sub {
  plan tests => 1;
  my $argument = "NC_000913.3.30k.fa";
  my $expected = "29442";
is ( retrieve_chromosome_size ($argument), $expected, 'retrieve chromosome size');
};

subtest 'write_chromosome_size_file' => sub {
  plan tests => 1;
  my @arguments = ( "assemblyHub/NC_000913.3/NC_000913.3.chrom.sizes", "NC_000913.3", "29930" );
  my $expected = "1";
  is ( write_chromosome_size_file (@arguments), $expected, 'write chromosome size file');
};

subtest 'convert_tracks' => sub {
  plan tests => 1;
  my @arguments = ( ".",
                    "assemblyHub/NC_000913.3",
                    "assemblyHub/NC_000913.3/NC_000913.3.chrom.sizes",
                    "log.txt" );
  my $expected = "1";
  is ( convert_tracks (@arguments), $expected, 'convert tracks');
};

# Check files
subtest 'Check files' => sub {
  plan tests => 12;
  # hub.txt
  file_md5sum_is ( "assemblyHub/hub.txt", "4b7e702e98019beb600a3493e55d37b6", "compare hub.txt files" );

  # genome.txt
  file_md5sum_is ( "assemblyHub/genome.txt", "06f05e4deb94d6ada2a5d42e02c7decb", "compare genome.txt files" );

  # fasta file
  file_md5sum_is ( "assemblyHub/NC_000913.3.fa", "608fa2e2ffbcf0104aab8d34c433b834", "compare fasta files" );

  # in $accession
  # .chrom.sizes
  file_md5sum_is ( "assemblyHub/NC_000913.3/NC_000913.3.chrom.sizes", "0782298d8959fa58881026197e878c30", "compare chrom.sizes files" );

  # .2bit file
  file_md5sum_is ( "assemblyHub/NC_000913.3/NC_000913.3.2bit", "674029770bfb21758a99fb8ed0e26ff7", "compare 2bit files" );

  # trackDb.txt
  file_md5sum_is ( "assemblyHub/NC_000913.3/trackDb.txt", "68b329da9893e34099c7d8ad5cb9c940", "compare trackDb.txt files" );

  # groups.txt
  file_md5sum_is ( "assemblyHub/NC_000913.3/groups.txt", "b653ebfb77f6b8eb29fb300d6e315f78", "compare groups.txt files" );

  # description.html
  file_md5sum_is ( "assemblyHub/NC_000913.3/description.html", "dedd73e9ee0a01169186aeb482120ad0", "compare description.html files" );

  # .bb files
  file_md5sum_is ( "assemblyHub/NC_000913.3/NC_000913.CDS.bb", "b30d0ca420966cbdedb24c093bbc762d", "compare bigbed files" );
  file_md5sum_is ( "assemblyHub/NC_000913.3/NC_000913.mobile_element.bb", "a89fbb5d2896abebaa3a2f7b5d8e0a24", "compare bigbed files" );
  file_md5sum_is ( "assemblyHub/NC_000913.3/NC_000913.ncRNA.bb", "94dbdf6f8d01c5ed947295fb24d0f0da", "compare bigbed files" );
  file_md5sum_is ( "assemblyHub/NC_000913.3/NC_000913.repeat_region.bb", "c90ef7a4b3534cb29fe48f34c22e512a", "compare bigbed files" );
};
