#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2016-12-06 17:12:00 mtw>

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Bio::ViennaNGS::Fasta;

my $fa_in = undef;
my $start = undef;
my $end = undef;
my @fids = ();


GetOptions('fasta|f=s'  => \$fa_in,
           'start|s=s'  => \$start,
           'end|e=s' => \$end);

unless (defined($fa_in)){die $!};
unless (defined($start)){die $!};
unless (defined($end)){die $!};

my $fastaO = Bio::ViennaNGS::Fasta->new(fa=>"$fa_in");

@fids = $fastaO->fastaids; # get all FASTA IDs

#print Dumper(\@fids);
#print ">>".ref(@fids)."<<\n";

foreach my $foo ($fids[0]){
  foreach my $id (@$foo){
    #print Dumper($id);

    my $seq = $fastaO->stranded_subsequence($id,
					  $start,
					  $end,
					  "+");
  print ">$id\n";
  print  join "\n", (unpack "(a70)*",$seq);
  print "\n";
  }
}

