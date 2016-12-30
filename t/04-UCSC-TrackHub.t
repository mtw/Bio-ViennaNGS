use strict;
use Bio::ViennaNGS::UCSC qw ( make_track_hub
                              retrieve_color
                              make_track
                              make_multi_bigwig_container_track
                              make_bigwig_container_track
                              retrieve_bigwig_tracks
                              retrieve_bigbed_url_tracks );
#use Test::Files;
use Test::File::Contents;
use File::Share ':all';
use Test::More tests => 12;
use Data::Dumper;

my $log = "TrackHub.log";

# Check if Bio::ViennaNGS::UCSC is successfully loaded.
require_ok ( 'Bio::ViennaNGS::UCSC');

# Make sure the module can do the functions.
can_ok ( 'Bio::ViennaNGS::UCSC', qw ( make_track_hub
                                      retrieve_color
                                      make_track
                                      make_multi_bigwig_container_track
                                      make_bigwig_container_track
                                      retrieve_bigwig_tracks
                                      retrieve_bigbed_url_tracks ));

# Check if functions work correctly
subtest 'retrieve_color' => sub {
  plan tests => 12;
  my @tests = 0..11;
  my @colours = ('133,154,0 133,154,0',
                '42,162,152 42,162,152',
                '38,140,210 38,140,210',
                '108,114,196 108,114,196',
                '211,55,130 211,55,130',
                '220,51,47 220,51,47',
                '203,76,22 203,76,22',
                '181,138,0 181,138,0',
                '0,44,54 0,44,54',
                '88,111,117 88,111,117',
                '133,154,0 133,154,0',
                '42,162,152 42,162,152');

  foreach my $number (@tests) {
      is ( retrieve_color($number), $colours[$number], "retrieve colour number $number");
  }
};

subtest 'make_tracks' => sub {
  plan tests => 2;
  my @arguments_pos = ( "hg19_highlyexpressed.pos.bb_bed",
                        "hg19_highlyexpressed.pos.bb_bed",
                        "http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.pos.bb",
                        "hg19_highlyexpressed",
                        "hg19_highlyexpressed",
                        "bigBed 12 .",
                        "off",
                        "Gene Id",
                        "name",
                        "133,154,0 133,154,0",
                        "pack",
                        "annotation",
                        "10" );
  my @arguments_neg = ( "hg19_highlyexpressed.neg.bb_bed",
                        "hg19_highlyexpressed.neg.bb_bed",
                        "http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.neg.bb",
                        "hg19_highlyexpressed",
                        "hg19_highlyexpressed",
                        "bigBed 12 .",
                        "off",
                        "Gene Id",
                        "name",
                        "42,162,152 42,162,152",
                        "pack",
                        "annotation",
                        "10" );
  my $expected_pos =
"#hg19_highlyexpressed.pos.bb_bed
track hg19_highlyexpressed.pos.bb_bed
bigDataUrl http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.pos.bb
shortLabel hg19_highlyexpressed
longLabel hg19_highlyexpressed
type bigBed 12 .
autoScale off
bedNameLabel Gene Id
searchIndex name
colorByStrand 133,154,0 133,154,0
visibility pack
group annotation
priority 10

";
  my $expected_neg =
"#hg19_highlyexpressed.neg.bb_bed
track hg19_highlyexpressed.neg.bb_bed
bigDataUrl http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.neg.bb
shortLabel hg19_highlyexpressed
longLabel hg19_highlyexpressed
type bigBed 12 .
autoScale off
bedNameLabel Gene Id
searchIndex name
colorByStrand 42,162,152 42,162,152
visibility pack
group annotation
priority 10

";
is ( make_track (@arguments_neg), $expected_neg, 'make neg. bigBed track');
is ( make_track (@arguments_pos), $expected_pos, 'make pos. bigBed track');
};

subtest 'make_multi_bigwig_container_track' => sub {
  plan tests => 1;
  my @arguments = ( "hg19_highlyexpressed_bw",
                    "hg19_highlyexpressed_bw",
                    "hg19_highlyexpressed_bw",
                    "hg19_highlyexpressed_bw",
                    "bigWig",
                    "on",
                    "full",
                    "1500" );
  my $expected =
"#hg19_highlyexpressed_bw
track hg19_highlyexpressed_bw
container multiWig
noInherit on
shortLabel hg19_highlyexpressed_bw
longLabel hg19_highlyexpressed_bw
type bigWig
configureable on
visibility full
aggregate transparentOverlay
showSubtrackColorOnUi on
autoScale on
windowingFunction maximum
priority 1500
alwaysZero on
yLineMark 0
yLineOnOff on
maxHeightPixels 125:125:11

";
is ( make_multi_bigwig_container_track (@arguments), $expected, 'make multi bigwig container track');
};

subtest 'make_pos_bigwig_container_track' => sub {
  plan tests => 1;
  my @arguments = ( "hg19_highlyexpressed_bw_pos",
                    "http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.pos.bw",
                    "hg19_highlyexpressed_bw_pos",
                    "hg19_highlyexpressed_bw_pos",
                    "bigWig",
                    "hg19_highlyexpressed_bw",
                    "133,154,0" );
  my $expected =
"track hg19_highlyexpressed_bw_pos
bigDataUrl http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.pos.bw
shortLabel hg19_highlyexpressed_bw_pos
longLabel hg19_highlyexpressed_bw_pos
type bigWig
parent hg19_highlyexpressed_bw
color 133,154,0

";
is ( make_bigwig_container_track (@arguments), $expected, 'make pos. bigwig container track');
};

subtest 'make_neg_bigwig_container_track' => sub {
  plan tests => 1;
  my @arguments = ( "hg19_highlyexpressed_bw_neg",
                    "http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.neg.bw",
                    "hg19_highlyexpressed_bw_neg",
                    "hg19_highlyexpressed_bw_neg",
                    "bigWig",
                    "hg19_highlyexpressed_bw",
                    "220,51,47" );
  my $expected =
"track hg19_highlyexpressed_bw_neg
bigDataUrl http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.neg.bw
shortLabel hg19_highlyexpressed_bw_neg
longLabel hg19_highlyexpressed_bw_neg
type bigWig
parent hg19_highlyexpressed_bw
color 220,51,47

";
is ( make_bigwig_container_track (@arguments), $expected, 'make neg bigwig container track');
};

subtest 'retrieve_bigwig_tracks' => sub {
  plan tests => 1;
  my @arguments = ( ".",
                    "http://www.tbi.univie.ac.at/~egg/hg19_trackHub/",
                    "trackHub",
                    "hg19",
                    "http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.pos.bw,http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.neg.bw" );
  my $expected =
"#hg19_highlyexpressed_bw
track hg19_highlyexpressed_bw
container multiWig
noInherit on
shortLabel hg19_highlyexpressed_bw
longLabel hg19_highlyexpressed_bw
type bigWig
configureable on
visibility full
aggregate transparentOverlay
showSubtrackColorOnUi on
autoScale on
windowingFunction maximum
priority 1500
alwaysZero on
yLineMark 0
yLineOnOff on
maxHeightPixels 125:125:11

track hg19_highlyexpressed_bw_pos
bigDataUrl http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.pos.bw
shortLabel hg19_highlyexpressed_bw_pos
longLabel hg19_highlyexpressed_bw_pos
type bigWig
parent hg19_highlyexpressed_bw
color 133,154,0

track hg19_highlyexpressed_bw_neg
bigDataUrl http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.neg.bw
shortLabel hg19_highlyexpressed_bw_neg
longLabel hg19_highlyexpressed_bw_neg
type bigWig
parent hg19_highlyexpressed_bw
color 220,51,47

";
is ( retrieve_bigwig_tracks (@arguments), $expected, 'retrieve bigwig tracks');
};

subtest 'retrieve_bigbed_url_tracks' => sub {
  plan tests => 1;
  my @arguments = ( ".",
                    "http://www.tbi.univie.ac.at/~egg/hg19_trackHub/",
                    "trackHub",
                    "http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.pos.bb#http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.neg.bb" );
  my $expected =
"#hg19_highlyexpressed.pos.bb_bed
track hg19_highlyexpressed.pos.bb_bed
bigDataUrl http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.pos.bb
shortLabel hg19_highlyexpressed
longLabel hg19_highlyexpressed
type bigBed 12 .
autoScale off
bedNameLabel Gene Id
searchIndex name
colorByStrand 133,154,0 133,154,0
visibility pack
group annotation
priority 10

#hg19_highlyexpressed.neg.bb_bed
track hg19_highlyexpressed.neg.bb_bed
bigDataUrl http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.neg.bb
shortLabel hg19_highlyexpressed
longLabel hg19_highlyexpressed
type bigBed 12 .
autoScale off
bedNameLabel Gene Id
searchIndex name
colorByStrand 42,162,152 42,162,152
visibility pack
group annotation
priority 10

";
is ( retrieve_bigbed_url_tracks (@arguments), $expected, 'retrieve bigbed url tracks');
};

# Make track hub
make_track_hub ( "hg19", ".", "http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_trackHub/", "http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.pos.bb#http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.neg.bb", "http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.pos.bw,http://nibiru.tbi.univie.ac.at/ViennaNGS/tutorial03/hg19_highlyexpressed.neg.bw", $log );

# Check files
# genome.txt
file_md5sum_is ( "trackHub/genome.txt", "336a16083ec54dade353cbc994122e9d", "compare hub.txt files");

# hub.txt
file_md5sum_is ( "trackHub/hub.txt", "fe4e2376c718c9edc5affa627734aae7", "compare genome.txt files" );

# trackDb.txt
file_md5sum_is ( "trackHub/hg19/trackDb.txt", "d35e774d83d9dc159d5d9679150c7ab7", "compare trackDb.txt files");
