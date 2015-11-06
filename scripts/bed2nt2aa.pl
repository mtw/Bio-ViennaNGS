#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2015-11-06 13:21:30 mtw>
#
# Provide nucleotide and amino acid sequences for BED files
#
# ***********************************************************************
# *  Copyright notice
# *
# *  Copyright 2015 Michael T. Wolfinger <michael@wolfinger.eu>
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
use Pod::Usage;
use Data::Dumper;
use Cwd;
use Path::Class;
use File::Slurp;
use Bio::ViennaNGS::Fasta;
use Bio::SeqUtils;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my ($fastaO,$bed,$line);
my $bed_in = undef;
my $fa_in  = undef;
my $translate = 0;  # translate nucleotide sequence

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("bed=s"      => \$bed_in,
					   "fa=s"       => \$fa_in,
					   "aa"         => sub{$translate=1},
					   "man"        => sub{pod2usage(-verbose => 2)},
					   "help|h"     => sub{pod2usage(1)}
					   );

unless(-f $bed_in){
  warn "Could not find input BED input file provided via ---bed option";
  pod2usage(-verbose => 0);
}
unless(-f $fa_in){
  warn "Could not find input FASTA input file provided via ---fa option";
  pod2usage(-verbose => 0);
}
my $cwd = getcwd();
unless ($bed_in =~ /^\// || $bed_in =~ /^\.\//){$bed_in = file($cwd,$bed_in);}
unless ($fa_in =~ /^\// || $fa_in =~ /^\.\//){$fa_in = file($cwd,$fa_in);}

$fastaO = Bio::ViennaNGS::Fasta->new(fa=>"$fa_in");

# parse BED6 file
$bed = read_file( $bed_in, array_ref => 1);
foreach $line (@$bed){
  chomp $line;
  my ($chr,$chromStart,$chromEnd,$name,$score,$strand) = split("\t",$line);
  my $seq = $fastaO->stranded_subsequence($chr,$chromStart+1,$chromEnd,$strand);
  print "$chr;$chromStart;$chromEnd;$strand;$seq;";
  if($translate == 1){
    my $seqobj = Bio::Seq->new (-seq => $seq,
				-id  => 'foobar',
				-accession_number => 'ABCD');
    my@prots = Bio::SeqUtils->translate_3frames($seqobj);
    foreach my $s (@prots){
      print $s->seq.";";
    }
    print "\n";
  }
}

__END__


=head1 NAME

bed2nt2aa.pl - Get nucleotde and amino acid sequences from BED intervals

=head1 SYNOPSIS

bed2nt2aa.pl [--bed I<FILE>] [--fa I<FILE>] [--aa] [options]

=head1 DESCRIPTION

Provide sequence information for BED intervals/files, optionally
translate nucleotide into amino acid sequences.

=head1 OPTIONS

=over

=item B<--bed>

Input file in BED6 format (mandatory)

=item B<--fa>

Input file in Fasta format (mandatory)

=item B<--aa>

Translate nucleotide into amino acid sequences, providing all three
possible frames

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut

