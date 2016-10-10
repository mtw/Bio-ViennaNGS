use strict;
use warnings;
use Path::Class;
use File::Share ':all';
use Test::More tests => 1;
use Bio::ViennaNGS::FeatureIO;
use Data::Dumper;

my $infile_path = dist_dir('Bio-ViennaNGS');
my $infile_bed12  = dist_file('Bio-ViennaNGS','data2/MIOS.bed12');

my @arg = (file => $infile_bed12,filetype => 'Bed12',instanceOf => 'Bed');
my $FIO = new_ok("Bio::ViennaNGS::FeatureIO" => \@arg);


