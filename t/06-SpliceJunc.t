# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl ViennaNGS.t'

#########################

# change 'tests => 1' to 'tests => last_test_to_print';

use strict;
use warnings;
#use IPC::Cmd qw(can_run);
use File::ShareDir::Install;
use Test::More tests => 24;
use Test::File::Contents;
use Data::Dumper;
use Path::Class;
use Cwd;
use Bio::ViennaNGS::SpliceJunc qw(bed6_ss_from_bed12 bed6_ss_from_rnaseq intersect_sj);
use Bio::ViennaNGS::Fasta;

my $cwd = getcwd();
my $data_dir_name = 'data-splicejunc-input';

my $data_dir = dir($cwd,$data_dir_name);

my $datasj_bed12_in = file($data_dir,'MIOS.bed12');
my $datasj_bed6_in = file($data_dir,'MIOS_sj.bed');
my $datasj_bed_ex = 'exist.SS.bed';
my $datasj_bed_nov = 'novel.SS.bed';
my $datasj_fasta_in = 'hg19_chromchecked.fa';
my %bed6_ss_from_bed12_results = (
  "chr7_7625266-7647110_uc010ktq.3.annotatedSS.bed6" => "c8e1d3564939f79606129f981b2cd303",
  "chr7_7622584-7647110_uc003srg.3.annotatedSS.bed6" => "b113664d78799c98dc8f2ed33f2a3ab9",
  "chr7_7607237-7613400_uc010ktp.1.annotatedSS.bed6" => "f74f051e613ac034a33cc6438a6fc90e",
  "chr7_7606615-7647110_uc003srf.3.annotatedSS.bed6" => "27b6132520902fa50ee1c8913a5080d1",
);
my %bed6_ss_from_rnaseq = (
  "chr7_7645702-7646625.mappedSS.bed6" => "953c09f527fb56816ddf0a1adc5ab8a1",
  "chr7_7636092-7645571.mappedSS.bed6" => "8188305a16ae31b2d02cfa21447904ab",
  "chr7_7634763-7635886.mappedSS.bed6" => "83ace65aa472e082677b9ce879f9ebcb",
  "chr7_7629194-7634609.mappedSS.bed6" => "ca927c5695b6795da56bbc0df567f890",
  "chr7_7628194-7634609.mappedSS.bed6" => "38d4dde16716bc0e9fca6d35fb8eeb4e",
  "chr7_7628194-7629034.mappedSS.bed6" => "81282645a2bc93883eab41ff20119c28",
  "chr7_7625436-7628127.mappedSS.bed6" => "b87979087db92ca8616d8f5dab8b15fa",
  "chr7_7623003-7625265.mappedSS.bed6" => "ef83be61cf34ff03015836b8b957a186",
  "chr7_7613827-7625265.mappedSS.bed6" => "68f07563b8e7baa1af74d462093df234",
  "chr7_7613827-7622747.mappedSS.bed6" => "703ec8fdd1e63e194f09784bc2021366",
  "chr7_7613400-7613727.mappedSS.bed6" => "0547196c3d27cf014e1dfc2e42c951e0",
  "chr7_7612224-7612336.mappedSS.bed6" => "ead51cb49e525be72ff0e3384384b559",
  "chr7_7607754-7612065.mappedSS.bed6" => "29363fe53392316c52a4446523f1691a",
  "chr7_7607319-7607655.mappedSS.bed6" => "94a11f9c38d61c42b38e5ade294f5bf0",
  "chr7_7607115-7607655.mappedSS.bed6" => "0025d7cd413cd4fac2f4878ee633c039",
  "chr7_7607002-7607655.mappedSS.bed6" => "0c4ed0b7f7af9c127e2bb40a13fcbd98",
  "chr7_7606703-7612065.mappedSS.bed6" => "16c3e486034157b786f3f974e1c3c96b",
  "chr7_7606703-7607655.mappedSS.bed6" => "b8cd47c18cbe451dcd31aaa407e415f8",
);
my %intersect_sj = (
  "novel.SS.bed" => "6239e270d534d7c6618f338de17e0112",
  "exist.SS.bed" => "85dd76401fddc563b0346a0808152b3e",
);

#########################

my $path_annot = './data-splicejunc-output/';
my $dest_ss = './data-splicejunc-output2/';
my $outdir = './data-splicejunc-output3/';
my $window = 10;
my $mincov = 0;
my $prefix = '';
my $max_intron_length = 7000;
my @result = ();
my $want_canonical = 0;
my $count = 1;
my $file = 0;

#########################

my $fastaO = Bio::ViennaNGS::Fasta->new(fa=>$datasj_fasta_in);

unless ($path_annot =~ /\/$/){$path_annot .= "/";}
unless (-d $path_annot){mkdir $path_annot or die $!;}

print STDERR "Testing routine bed6_ss_from_bed12...\n";
bed6_ss_from_bed12($datasj_bed12_in,$path_annot,$window,$want_canonical,$fastaO);

foreach my $el (keys %bed6_ss_from_bed12_results) {
  file_md5sum_is $path_annot.$el, $bed6_ss_from_bed12_results{$el}, "output $count is the same";
  $count++;
}


print STDERR "\nTesting routine bed6_ss_from_rnaseq...\n";
unless ($dest_ss =~ /\/$/){$dest_ss .= "/";}
unless (-d $dest_ss){mkdir $dest_ss or die $!;}

bed6_ss_from_rnaseq($datasj_bed6_in,$dest_ss,$window,$mincov,$want_canonical,$fastaO);


foreach my $el (keys %bed6_ss_from_rnaseq) {
  file_md5sum_is $dest_ss.$el, $bed6_ss_from_rnaseq{$el}, "output $count is the same";
  $count++;
}


print STDERR "\nTesting routine intersect_sj...\n";
unless ($outdir =~ /\/$/){$outdir .= "/";}
unless (-d $outdir){mkdir $outdir or die $!;}

@result = intersect_sj($path_annot,$dest_ss,$outdir,$prefix,$window,$max_intron_length);
my ($exist,$novel) = @result;

file_md5sum_is $outdir.$datasj_bed_ex, $intersect_sj{$datasj_bed_ex}, "output 23 is the same";
file_md5sum_is $outdir.$datasj_bed_nov, $intersect_sj{$datasj_bed_nov}, "output 24 is the same";




#my @data1_ids = $f->fastaids;
#ok (scalar @data1_ids, 1);
#my $ps = $f->primaryseq;


# Insert your test code below, the Test::More module is use()ed here so read
# its man page ( perldoc Test::More ) for help writing this test script.



done_testing;
