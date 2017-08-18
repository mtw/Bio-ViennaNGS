#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2017-08-17 18:11:05 mtw>
#
# Find sequence motifs in (multi)-Fasta files
#
# usage: fasta_regex.pl --motif <STRING> --fa <FASTAFILE>
#
# ***********************************************************************
# *  Copyright notice
# *
# *  Copyright 2017 Michael T. Wolfinger <michael@wolfinger.eu>
# *  All rights reserved
# *
# *  This program is free software: you can redistribute it and/or modify
# *  it under the terms of the GNU General Public License as published by
# *  the Free Software Foundation, either version 3 of the License, or
# *  (at your option) any later version.
# *
# *  This program is distributed in the hope that it will be useful,
# *  but WITHOUT ANY WARRANTY; without even the implied warranty of
# *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# *  GNU General Public License for more details.
# *
# *  You should have received a copy of the GNU General Public License
# *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
# *
# *  This copyright notice MUST APPEAR in all copies of the script!
# ***********************************************************************

use strict;
use warnings;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Data::Dumper;
use Pod::Usage;
use Cwd;
use Path::Class;
use Bio::ViennaNGS::Fasta;
use Carp;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my ($fastaO,$cwd);
my $motif  = undef;
my $fa_in  = undef;
my $strand = "+";
my $wantbed = 0;
my $name   = "pattern";
my @fids = ();


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^ Main ^^^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("fa|f=s"      => \$fa_in,
					   "motif|m=s"   => \$motif,
					   "name|n=s"    => \$name,
					   "bed"         => sub{$wantbed = 1},
					   "man"         => sub{pod2usage(-verbose => 2)},
					   "help|h"      => sub{pod2usage(1)}
					  );

unless (-f $fa_in){
  warn "Could not find input FASTA input file provided via ---fa option";
  pod2usage(-verbose => 0);
}
$cwd = getcwd();
unless ($fa_in =~ /^\// || $fa_in =~ /^\.\//){$fa_in = file($cwd,$fa_in);}

unless (defined $motif) {
  warn "No search pattern found. Please provided a mtoif as regular expression.";
  pod2usage(-verbose => 0);
}
$strand eq '-' ? ($strand = '-') : ($strand = '+') ;

$fastaO = Bio::ViennaNGS::Fasta->new(fasta=>"$fa_in");
my $ps = $fastaO->primaryseqH;
foreach my $id (keys %$ps){
  my $cnt = 0;
  my $start = $$ps{$id}->{start};
  my $end = $$ps{$id}->{stop};
 # print ">FASTA Id $id ($start-$end)\n";
  my $seq = $fastaO->stranded_subsequence($id,$start, $end,$strand);
 # print "$seq\n";
  while ($seq =~ m/($motif)/g){ # we need those parentheses here !!!
    my $p = pos($seq)-length($1);
    my $q = $p+length($1);
    my $n = $name."_$p-$q";
    if ($wantbed == 1){
      print join "\t", ($id,$p,$q,$n,"0",$strand);
      print "\n";
    }
    $cnt++;
  }
  #print "$cnt motifs found\n";
}


#print Dumper($ps);

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^^ Subroutines ^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#




__END__

=head1 NAME

fasta_regex.pl - Find sequence motifs in a (multi) Fasta file

=head1 SYNOPSIS

fasta_regex.pl [--motif I<REGEX>] [--fa I<FILE>] [--name I<STRING>] [--bed]
[options]

=head1 DESCRIPTION

This script extracts sequence motifs, provided as regular expressionm
from Fasta files.

The tool returns B<all motifs> matching the search criteria as a
(mutli-)Fasta file to STDOUT. This means if the motif is found more
than once within an annotated feature, all matches will be
reported. Optionally, all loci matching the regular expresion pattern
are reorted as BED6 intervals.

=head1 OPTIONS

=over

=item B<fa|f>

Reference genome in Fasta format

=item B<motif|m>

The motif to search for as regular expression.

=item B<bed>

Switch for printing BED6 intervals of matching motif loci.

=item B<name|n>

Optional name for BED6 intervals.

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut
