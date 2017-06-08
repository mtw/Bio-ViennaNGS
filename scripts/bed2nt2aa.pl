#!/usr/bin/env perl
# -*-CPerl-*-
# Last changed Time-stamp: <2017-06-07 13:43:06 michl>
#
# Extract nucleotide and amino acid sequence data from Fasta file,
# provided by a BED file
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
use Pod::Usage;
use Data::Dumper;
use Cwd;
use Carp;
use Path::Class;
use Bio::ViennaNGS::Fasta;
use Bio::ViennaNGS::FeatureIO;
use diagnostics;

use Bio::SeqUtils;

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
#^^^^^^^^^^ Variables ^^^^^^^^^^^#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#

my ($fastaO,$bed,$line,$cwd,$f,$seq,$theid,$intervals,$out);
my $bed_in = undef;
my $fa_in  = undef;
my $stdout = 0;
my $translate = 0;  # translate nucleotide sequence

Getopt::Long::config('no_ignore_case');
pod2usage(-verbose => 1) unless GetOptions("bed=s"      => \$bed_in,
					   "fa=s"       => \$fa_in,
					   "a|aa"       => sub{$translate=1},
					   "s|stdout"   => sub{$stdout=1},
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
unless ($fa_in  =~ /^\// || $fa_in  =~ /^\.\//){$fa_in  = file($cwd,$fa_in);}

$fastaO    = Bio::ViennaNGS::Fasta->new(fasta=>"$fa_in");
$intervals = Bio::ViennaNGS::FeatureIO->new(
					    file => "$bed_in",
					    filetype => 'Bed6',
					    instanceOf => 'Feature',
					   );
foreach  $f (@{$intervals->data}){
  confess "ERROR: ID ".$f->chromosome." not found in Fasta file $fa_in"
    unless ($fastaO->has_sequid($f->chromosome));
  $seq = $fastaO->stranded_subsequence($f->chromosome,
				       $f->start+1,
				       $f->end,
				       $f->strand);
  #my $fn = $f->chromosome.".".$f->name.".".eval($f->start+1)."-".$f->end.".fa";
  my $fn = "sequence.fa";
  my $newid = $f->chromosome."|".$f->name."|".eval($f->start+1)."-".$f->end."|";
  $stdout == 0 ? open($out, ">", $fn) : open($out, ">&STDOUT");
  print $out ">$newid\n";
  print $out join "\n", (unpack "(a70)*",$seq);
  print $out "\n";
  if($translate == 1){
    my $seqobj = Bio::Seq->new (-seq => $seq,
				-id  => 'foobar',
				-accession_number => 'ABCD');
    my @prots = Bio::SeqUtils->translate_3frames($seqobj);
    my $i=1;
    foreach my $s (@prots){
      print $out ">$newid|AA frame $i|\n"; $i++;
      print $out join "\n", (unpack "(a70)*",$s->seq);
      print $out "\n";
    }
  }
  close $out;
}

__END__


=head1 NAME

bed2nt2aa.pl - Get nucleotde and amino acid sequences from BED intervals

=head1 SYNOPSIS

bed2nt2aa.pl [--bed I<FILE>] [--fa I<FILE>] [--aa] [options]

=head1 DESCRIPTION

Extract sequence intervals from Fasta files, provided as BED6
records. Optionally translate nucleotide into amino acid
sequences. Output is written to a 'sequence.fa' output file. This can
be overriden with the C<-s> switch, which triggers STDOUT output.

=head1 OPTIONS

=over

=item B<--bed>

Input file in BED6 format (mandatory)

=item B<--fa>

Input file in Fasta format (mandatory)

=item B<-a|--aa>

Translate nucleotide into amino acid sequences, providing all three
possible frames

=item B<-s|--stdout>

Write output to STDOUT instead of an output file 'sequence.fa'

=item B<--help -h>

Print short help

=item B<--man>

Prints the manual page and exits

=back

=head1 AUTHOR

Michael T. Wolfinger E<lt>michael@wolfinger.euE<gt>

=cut

