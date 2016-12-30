#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2016-11-02 19:36:08 mtw>
#
# Provide nucleotide and amino acid sequences for BED files
#
# ***********************************************************************
# *  Copyright notice
# *
# *  Copyright 2016 Michael T. Wolfinger <michael@wolfinger.eu>
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
use Bio::ViennaNGS::Fasta;
use Bio::ViennaNGS::FeatureIO;

use Bio::SeqUtils;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my ($fastaO,$bed,$line,$cwd,$f,$seq,$theid, $intervals);
my @fids = ();
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
$cwd = getcwd();
unless ($bed_in =~ /^\// || $bed_in =~ /^\.\//){$bed_in = file($cwd,$bed_in);}
unless ($fa_in =~ /^\// || $fa_in =~ /^\.\//){$fa_in = file($cwd,$fa_in);}

$fastaO = Bio::ViennaNGS::Fasta->new(fa=>"$fa_in");
@fids = $fastaO->fastaids; # get all FASTA IDs
$theid = $fids[0][0];
$intervals = Bio::ViennaNGS::FeatureIO->new(
					    file => "$bed_in",
					    filetype => 'Bed6',
					    instanceOf => 'Feature',
					   );
foreach  $f (@{$intervals->data}){
  $seq = $fastaO->stranded_subsequence($f->chromosome,
				       $f->start+1,
				       $f->end,
				       $f->strand);
  my $fn = $f->chromosome.".".$f->name.".".eval($f->start+1)."-".$f->end.".fa";
  my $newid = $theid."|".$f->name."|".eval($f->start+1)."-".$f->end."|";
  open (FH, ">", $fn);
  print FH ">$newid\n";
  print FH join "\n", (unpack "(a70)*",$seq);
  print FH "\n";
  if($translate == 1){
    my $seqobj = Bio::Seq->new (-seq => $seq,
				-id  => 'foobar',
				-accession_number => 'ABCD');
    my @prots = Bio::SeqUtils->translate_3frames($seqobj);
    my $i=1;
    foreach my $s (@prots){
      print FH "\n";
      print FH ">$newid|AA frame $i|\n"; $i++;
      print FH join "\n", (unpack "(a70)*",$s->seq);
      print FH "\n";
    }
  }
  close FH;
}

__END__


=head1 NAME

bed2nt2aa.pl - Get nucleotde and amino acid sequences from BED intervals

=head1 SYNOPSIS

bed2nt2aa.pl [--bed I<FILE>] [--fa I<FILE>] [--aa] [options]

=head1 DESCRIPTION

Generate sequence information for BED6 intervals as Fasta files, optionally
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

